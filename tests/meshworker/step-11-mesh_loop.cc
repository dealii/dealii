// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// step-11 but rewritten using mesh_loop()


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// Just this one is new: it declares a class
// DynamicSparsityPattern, which we will use and explain
// further down below.
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// We will make use of the std::find algorithm of the C++ standard library, so
// we have to include the following file for its declaration:
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

// The last step is as in all previous programs:
namespace Step11
{
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
    LaplaceProblem(const unsigned int mapping_degree);
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_and_solve();
    void
    solve();

    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;
    MappingQ<dim>      mapping;

    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    AffineConstraints<double> constraints;

    Vector<double> solution;
    Vector<double> system_rhs;

    TableHandler output_table;
  };



  // Construct such an object, by initializing the variables. Here, we use
  // linear finite elements (the argument to the <code>fe</code> variable
  // denotes the polynomial degree), and mappings of given order. Print to
  // screen what we are about to do.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem(const unsigned int mapping_degree)
    : fe(1)
    , dof_handler(triangulation)
    , mapping(mapping_degree)
  {
    deallog << "Using mapping with degree " << mapping_degree << ":\n"
            << "============================" << std::endl;
  }



  // The first task is to set up the variables for this problem. This includes
  // generating a valid <code>DoFHandler</code> object, as well as the
  // sparsity patterns for the matrix, and the object representing the
  // constraints that the mean value of the degrees of freedom on the boundary
  // be zero.
  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
    DoFTools::extract_boundary_dofs(dof_handler,
                                    ComponentMask(),
                                    boundary_dofs);

    const unsigned int first_boundary_dof = std::distance(
      boundary_dofs.begin(),
      std::find(boundary_dofs.begin(), boundary_dofs.end(), true));

    // Then generate a constraints object with just this one constraint. First
    // clear all previous content (which might reside there from the previous
    // computation on a once coarser grid), then add this one line
    // constraining the <code>first_boundary_dof</code> to the sum of other
    // boundary DoFs each with weight -1. Finally, close the constraints
    // object, i.e. do some internal bookkeeping on it for faster processing
    // of what is to come later:
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.add_line(first_boundary_dof);
    for (unsigned int i = first_boundary_dof + 1; i < dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true)
        constraints.add_entry(first_boundary_dof, i, -1);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);


    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }



  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim> &      mapping,
                const FiniteElement<dim> &fe,
                const unsigned int        quadrature_degree)
      : fe_values(mapping,
                  fe,
                  QGauss<dim>(quadrature_degree),
                  update_values | update_gradients | update_quadrature_points |
                    update_JxW_values)
      , fe_face_values(mapping,
                       fe,
                       QGauss<dim - 1>(quadrature_degree),
                       update_values | update_quadrature_points |
                         update_JxW_values | update_normal_vectors)
    {}

    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  update_values | update_gradients | update_quadrature_points |
                    update_JxW_values)
      , fe_face_values(scratch_data.fe_face_values.get_mapping(),
                       scratch_data.fe_face_values.get_fe(),
                       scratch_data.fe_face_values.get_quadrature(),
                       update_values | update_quadrature_points |
                         update_JxW_values | update_normal_vectors)
    {}

    FEValues<dim>     fe_values;
    FEFaceValues<dim> fe_face_values;
  };
  struct CopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };

  // The next function then assembles the linear system of equations, solves
  // it, and evaluates the solution.
  template <int dim>
  void
  LaplaceProblem<dim>::assemble_and_solve()
  {
    typedef decltype(dof_handler.begin_active()) Iterator;

    auto cell_worker = [](const Iterator &  cell,
                          ScratchData<dim> &scratch_data,
                          CopyData &        copy_data) {
      const unsigned int dofs_per_cell =
        scratch_data.fe_values.get_fe().dofs_per_cell;
      const unsigned int n_q_points =
        scratch_data.fe_values.get_quadrature().size();

      copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      copy_data.cell_rhs.reinit(dofs_per_cell);

      copy_data.local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(copy_data.local_dof_indices);

      scratch_data.fe_values.reinit(cell);

      std::vector<double>   rhs_values(n_q_points);
      ConstantFunction<dim> right_hand_side(-2.0);
      right_hand_side.value_list(scratch_data.fe_values.get_quadrature_points(),
                                 rhs_values);

      // ... and assemble the local contributions to the system matrix and
      // right hand side as also discussed above:
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=
                (scratch_data.fe_values.shape_grad(i, q_point) *
                 scratch_data.fe_values.shape_grad(j, q_point) *
                 scratch_data.fe_values.JxW(q_point));

            copy_data.cell_rhs(i) +=
              (scratch_data.fe_values.shape_value(i, q_point) *
               rhs_values[q_point] * scratch_data.fe_values.JxW(q_point));
          }
    };

    auto boundary_worker = [](const Iterator &    cell,
                              const unsigned int &face_no,
                              ScratchData<dim> &  scratch_data,
                              CopyData &          copy_data) {
      const unsigned int dofs_per_cell =
        scratch_data.fe_values.get_fe().dofs_per_cell;
      const unsigned int n_face_q_points =
        scratch_data.fe_face_values.get_quadrature().size();

      std::vector<double>   face_boundary_values(n_face_q_points);
      ConstantFunction<dim> boundary_values(1.0);

      scratch_data.fe_face_values.reinit(cell, face_no);
      boundary_values.value_list(
        scratch_data.fe_face_values.get_quadrature_points(),
        face_boundary_values);

      for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          copy_data.cell_rhs(i) -=
            (face_boundary_values[q_point] *
             scratch_data.fe_face_values.shape_value(i, q_point) *
             scratch_data.fe_face_values.JxW(q_point));
    };

    auto copier = [&](const CopyData &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             system_rhs);
    };

    const unsigned int gauss_degree =
      std::max(static_cast<unsigned int>(
                 std::ceil(1. * (mapping.get_degree() + 1) / 2)),
               2U);

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          ScratchData<dim>(mapping, fe, gauss_degree),
                          CopyData(),
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces,
                          boundary_worker);


    solve();


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
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      ZeroFunction<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(gauss_degree + 1),
                                      VectorTools::H1_seminorm);
    // Then, the function just called returns its results as a vector of
    // values each of which denotes the norm on one cell. To get the global
    // norm, we do the following:
    const double norm =
      VectorTools::compute_global_error(triangulation,
                                        norm_per_cell,
                                        VectorTools::H1_seminorm);

    // Last task -- generate output:
    output_table.add_value("cells", triangulation.n_active_cells());
    output_table.add_value("|u|_1", norm);
    output_table.add_value("error",
                           std::fabs(norm - std::sqrt(3.14159265358 / 2)));
  }



  // The following function solving the linear system of equations is copied
  // from step-5 and is explained there in some detail:
  template <int dim>
  void
  LaplaceProblem<dim>::solve()
  {
    SolverControl solver_control(1000, 1e-12, false, false);
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }



  // Finally the main function controlling the different steps to be
  // performed.
  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    GridGenerator::hyper_ball(triangulation);
    static const SphericalManifold<dim> boundary;
    triangulation.set_all_manifold_ids_on_boundary(0);
    triangulation.set_manifold(0, boundary);

    for (unsigned int cycle = 0; cycle < 5; ++cycle)
      {
        if (cycle > 0)
          triangulation.refine_global(1);

        setup_system();
        assemble_and_solve();
      };

    // After all the data is generated, write a table of results to the
    // screen:
    output_table.set_precision("|u|_1", 6);
    output_table.set_precision("error", 6);
    output_table.write_text(deallog.get_file_stream());
    deallog << std::endl;
  }
} // namespace Step11

int
main()
{
  initlog();

  for (unsigned int mapping_degree = 1; mapping_degree <= 3; ++mapping_degree)
    Step11::LaplaceProblem<2>(mapping_degree).run();

  return 0;
}
