/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2013 by the deal.II authors
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
 * Author: Markus Buerg, University of Karlsruhe, 2010
 */


// @sect3{Include files}

// The include files are already known. The one critical for the current
// program is the one that contains the ConstraintMatrix in the
// <code>lac/</code> directory:
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>


namespace Step45
{
  using namespace dealii;

  // @sect3{The <code>LaplaceProblem</code> class}

  // The class <code>LaplaceProblem</code> is the main class of this
  // problem. As mentioned in the introduction, it is fashioned after the
  // corresponding class in step-3. Correspondingly, the documentation from
  // that tutorial program applies here as well. The only new member variable
  // is the <code>constraints</code> variables that will hold the constraints
  // from the periodic boundary condition. We will initialize it in the
  // <code>make_periodicity_constraints()</code> function which we call from
  // <code>make_grid_and_dofs()</code>.
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    void run ();

  private:
    Triangulation<2> triangulation;

    FE_Q<2> fe;
    DoFHandler<2> dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> system_rhs;
    Vector<double> solution;

    void assemble_system ();
    void output_results ();
    void make_grid_and_dofs ();
    void make_periodicity_constraints ();
    void solve ();
  };


  // @sect3{The <code>RightHandSide</code> class}

  // The following implements the right hand side function discussed in the
  // introduction. Its implementation is obvious given what has been shown in
  // step-4 before:
  class RightHandSide: public Function<2>
  {
  public:
    RightHandSide ();

    virtual double value (const Point<2> &p,
                          const unsigned int component = 0) const;
  };


  RightHandSide::RightHandSide ()
    :
    Function<2> ()
  {}


  double
  RightHandSide::value (const Point<2> &p,
                        const unsigned int) const
  {
    return (std::cos (2 * numbers::PI * p(0)) *
            std::exp (- 2 * p(0)) *
            std::cos (2 * numbers::PI * p(1)) *
            std::exp (- 2 * p(1)));
  }

  // @sect3{Implementation of the <code>LaplaceProblem</code> class}

  // The first part of implementing the main class is the constructor. It is
  // unchanged from step-3 and step-4:
  LaplaceProblem::LaplaceProblem ()
    :
    fe (1),
    dof_handler (triangulation)
  {}


  // @sect4{LaplaceProblem::make_grid_and_dofs}

  // The following is the first function to be called in
  // <code>run()</code>. It sets up the mesh and degrees of freedom.
  //
  // We start by creating the usual square mesh and changing the boundary
  // indicator on the parts of the boundary where we have Dirichlet boundary
  // conditions (top and bottom, i.e. faces two and three of the reference
  // cell as defined by GeometryInfo), so that we can distinguish between the
  // parts of the boundary where periodic and where Dirichlet boundary
  // conditions hold. We then refine the mesh a fixed number of times, with
  // child faces inheriting the boundary indicators previously set on the
  // coarse mesh from their parents.
  void LaplaceProblem::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation);
    triangulation.begin_active ()->face (2)->set_boundary_indicator (1);
    triangulation.begin_active ()->face (3)->set_boundary_indicator (1);
    triangulation.refine_global (5);

    // The next step is to distribute the degrees of freedom and produce a
    // little bit of graphical output:
    dof_handler.distribute_dofs (fe);
    std::cout << "Number of active cells: "
              << triangulation.n_active_cells ()
              << std::endl
              << "Degrees of freedom: " << dof_handler.n_dofs ()
              << std::endl;

    // Now it is the time for the constraints that come from the periodicity
    // constraints. We do this in the following, separate function, after
    // clearing any possible prior content from the constraints object:
    constraints.clear ();
    make_periodicity_constraints ();

    // We also incorporate the homogeneous Dirichlet boundary conditions on
    // the upper and lower parts of the boundary (i.e. the ones with boundary
    // indicator 1) and close the <code>ConstraintMatrix</code> object:
    VectorTools::interpolate_boundary_values (dof_handler, 1,
                                              ZeroFunction<2> (),
                                              constraints);
    constraints.close ();

    // Then we create the sparsity pattern and the system matrix and
    // initialize the solution and right-hand side vectors. This is again as
    // in step-3 or step-6, for example:
    CompressedSparsityPattern c_sparsity_pattern (dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler,
                                     c_sparsity_pattern,
                                     constraints,
                                     false);
    c_sparsity_pattern.compress ();
    sparsity_pattern.copy_from (c_sparsity_pattern);

    system_matrix.reinit (sparsity_pattern);
    system_rhs.reinit (dof_handler.n_dofs());
    solution.reinit (dof_handler.n_dofs());
  }



  // @sect4{LaplaceProblem::make_periodicity_constraints}

  // This is the function that provides the new material of this tutorial
  // program. The general outline of the algorithm is as follows: we first
  // loop over all the degrees of freedom on the right boundary and record
  // their $y$-locations in a map together with their global indices. Then we
  // go along the left boundary, find matching $y$-locations for each degree
  // of freedom, and then add constraints that identify these matched degrees
  // of freedom.
  //
  // In this function, we make use of the fact that we have a scalar element
  // (i.e. the only valid vector component that can be passed to
  // DoFAccessor::vertex_dof_index is zero) and that we have a $Q_1$ element
  // for which all degrees of freedom live in the vertices of the
  // cell. Furthermore, we have assumed that we are in 2d and that meshes were
  // not refined adaptively &mdash; the latter assumption would imply that
  // there may be vertices that aren't matched one-to-one and for which we
  // won't be able to compute constraints this easily. We will discuss in the
  // "outlook" part of the results section below other strategies to write the
  // current function that can work in cases like this as well.
  void LaplaceProblem::make_periodicity_constraints ()
  {
    // To start with the actual implementation, we loop over all active cells
    // and check whether the cell is located at the right boundary (i.e. face
    // 1 &mdash; the one at the right end of the cell &mdash; is at the
    // boundary). If that is so, then we use that for the currently used
    // finite element, each degree of freedom of the face is located on one
    // vertex, and store their $y$-coordinate along with the global number of
    // this degree of freedom in the following map:
    std::map<unsigned int, double> dof_locations;

    for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active ();
         cell != dof_handler.end (); ++cell)
      if (cell->at_boundary ()
          &&
          cell->face(1)->at_boundary ())
        {
          dof_locations[cell->face(1)->vertex_dof_index(0, 0)]
            = cell->face(1)->vertex(0)[1];
          dof_locations[cell->face(1)->vertex_dof_index(1, 0)]
            = cell->face(1)->vertex(1)[1];
        }
    // Note that in the above block, we add vertices zero and one of the
    // affected face to the map. This means that we will add each vertex
    // twice, once from each of the two adjacent cells (unless the vertex is a
    // corner of the domain). Since the coordinates of the vertex are the same
    // both times of course, there is no harm: we replace one value in the map
    // with itself the second time we visit an entry.
    //
    // The same will be true below where we add the same constraint twice to
    // the ConstraintMatrix &mdash; again, we will overwrite the constraint
    // with itself, and no harm is done.

    // Now we have to find the corresponding degrees of freedom on the left
    // part of the boundary. Therefore we loop over all cells again and choose
    // the ones where face 0 is at the boundary:
    for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active ();
         cell != dof_handler.end (); ++cell)
      if (cell->at_boundary ()
          &&
          cell->face (0)->at_boundary ())
        {
          // Every degree of freedom on this face needs to have a
          // corresponding one on the right side of the face, and our goal is
          // to add a constraint for the one on the left in terms of the one
          // on the right. To this end we first add a new line to the
          // constraint matrix for this one degree of freedom. Then we
          // identify it with the corresponding degree of freedom on the right
          // part of the boundary by constraining the degree of freedom on the
          // left with the one on the right times a weight of 1.0.
          //
          // Consequently, we loop over the two vertices of each face we find
          // and then loop over all the $y$-locations we've previously
          // recorded to find which degree of freedom on the right boundary
          // corresponds to the one we currently look at. Note that we have
          // entered these into a map, and when looping over the iterators
          // <code>p</code> of this map, <code>p-@>first</code> corresponds to
          // the "key" of an entry (the global number of the degree of
          // freedom), whereas <code>p-@>second</code> is the "value" (the
          // $y$-location we have entered above).
          //
          // We are quite sure here that we should be finding such a
          // corresponding degree of freedom. However, sometimes stuff happens
          // and so the bottom of the block contains an assertion that our
          // assumption was indeed correct and that a vertex was found.
          for (unsigned int face_vertex = 0; face_vertex<2; ++face_vertex)
            {
              constraints.add_line (cell->face(0)->vertex_dof_index (face_vertex, 0));

              std::map<unsigned int, double>::const_iterator p = dof_locations.begin();
              for (; p != dof_locations.end(); ++p)
                if (std::fabs(p->second - cell->face(0)->vertex(face_vertex)[1]) < 1e-8)
                  {
                    constraints.add_entry (cell->face(0)->vertex_dof_index (face_vertex, 0),
                                           p->first, 1.0);
                    break;
                  }
              Assert (p != dof_locations.end(),
                      ExcMessage ("No corresponding degree of freedom was found!"));
            }
        }
  }



  // @sect4{LaplaceProblem::assemble_system}

  // Assembling the system matrix and the right-hand side vector is done as in
  // other tutorials before.
  //
  // The only difference here is that we don't copy elements from local
  // contributions into the global matrix and later fix up constrained degrees
  // of freedom, but that we let the ConstraintMatrix do this job in one swoop
  // for us using the ConstraintMatrix::distribute_local_to_global
  // function(). This was previously already demonstrated in step-16, step-22,
  // for example, along with a discussion in the introduction of step-27.
  void LaplaceProblem::assemble_system ()
  {
    QGauss<2>  quadrature_formula(2);
    FEValues<2> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const RightHandSide right_hand_side;

    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                        endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;
        cell_rhs = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                     fe_values.shape_grad (j, q_point) *
                                     fe_values.JxW (q_point));

              cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                              right_hand_side.value (fe_values.quadrature_point (q_point)) *
                              fe_values.JxW (q_point));
            }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);
      }
  }


  // @sect4{LaplaceProblem::solve}

  // To solve the linear system of equations $Au=b$ we use the CG solver with
  // an SSOR-preconditioner. This is, again, copied almost verbatim from
  // step-6. As in step-6, we need to make sure that constrained degrees of
  // freedom get their correct values after solving by calling the
  // ConstraintMatrix::distribute function:
  void LaplaceProblem::solve ()
  {
    SolverControl solver_control (dof_handler.n_dofs (), 1e-12);
    PreconditionSSOR<SparseMatrix<double> > precondition;

    precondition.initialize (system_matrix);

    SolverCG<> cg (solver_control);

    cg.solve (system_matrix, solution, system_rhs, precondition);
    constraints.distribute (solution);
  }


  // @sect4{LaplaceProblem::output_results}

  // This is another function copied from previous tutorial programs. It
  // generates graphical output in VTK format:
  void LaplaceProblem::output_results ()
  {
    DataOut<2> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "u");
    data_out.build_patches ();

    std::ofstream output ("solution.vtk");

    data_out.write_vtk (output);
  }



  // @sect4{LaplaceProblem::run}

  // And another function copied from previous programs:
  void LaplaceProblem::run ()
  {
    make_grid_and_dofs();
    assemble_system ();
    solve ();
    output_results ();
  }
}

// @sect3{The <code>main</code> function}

// And at the end we have the main function as usual, this time copied from
// step-6:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step45;

      deallog.depth_console (0);

      LaplaceProblem laplace_problem;
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
