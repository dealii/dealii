// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/**
 * @defgroup assembly Assembly of finite element operators
 * @ingroup Integrators
 *
 * Finite element operators are built from contributions on cells and, for
 * some formulations, on faces. deal.II offers several ways to organize this
 * computation. The choice concerns both the traversal of the mesh and whether
 * the global matrix is stored explicitly.
 *
 * <h3>Explicit cell loops</h3>
 *
 * In step-3, an explicit loop visits each active cell. FEValues evaluates
 * shape functions and quadrature data, nested loops compute a local matrix
 * and right hand side, and their entries are added to the global objects.
 * This exposes all parts of assembly and is a useful starting point for a
 * new formulation. Applications with hanging nodes or other constraints can
 * use AffineConstraints::distribute_local_to_global() to transfer local
 * contributions, as shown in step-6.
 *
 * For example, the cell loop for a scalar Laplace matrix looks as follows.
 * Here @p fe is the finite element, @p constraints is closed, and
 * @p system_matrix has the required size and sparsity pattern. FEValues
 * requests gradients and quadrature weights; the local objects are reused
 * between cells.
 * @code
 * FEValues<dim> fe_values(fe, QGauss<dim>(fe.degree + 1),
 *                         update_gradients | update_JxW_values);
 * const unsigned int n = fe.n_dofs_per_cell();
 * FullMatrix<double> cell_matrix(n, n);
 * std::vector<types::global_dof_index> local_dof_indices(n);
 * for (const auto &cell : dof_handler.active_cell_iterators())
 *   {
 *     fe_values.reinit(cell);
 *     cell_matrix = 0;
 *     for (const unsigned int q : fe_values.quadrature_point_indices())
 *       for (unsigned int i = 0; i < n; ++i)
 *         for (unsigned int j = 0; j < n; ++j)
 *           cell_matrix(i, j) += fe_values.shape_grad(i, q) *
 *                                fe_values.shape_grad(j, q) * fe_values.JxW(q);
 *     cell->get_dof_indices(local_dof_indices);
 *     constraints.distribute_local_to_global(cell_matrix, local_dof_indices,
 *                                            system_matrix);
 *   }
 * @endcode
 *
 * <h3>WorkStream</h3>
 *
 * WorkStream::run() separates the computation of local contributions from
 * their insertion into global objects. A worker uses scratch data to compute
 * each cell's contribution and stores the result in copy data. A copier then
 * inserts that result into the global matrix and right hand side. This
 * separation allows local work to run on multiple threads. For the
 * iterator-range overload, insertion into shared objects is serialized.
 * The colored overload can also run copiers concurrently: its coloring must
 * ensure that copiers within one color do not write to the same memory.
 * The local finite element calculations can remain the same as in an explicit
 * cell loop. See step-9 and the
 * @ref threads "shared-memory parallelism documentation" for examples and
 * details of the WorkStream model.
 *
 * The iterator-range overload replaces the explicit cell loop by a worker
 * and a copier. In this sketch, @p ScratchData owns an FEValues object and
 * @p CopyData owns the local matrix and DoF indices. The application function
 * @p assemble_cell resets and fills the copy data using the local calculation
 * above; it only writes to its scratch and copy arguments. Both data types
 * must satisfy WorkStream's copying requirements (see step-9).
 * @code
 * using Iterator = typename DoFHandler<dim>::active_cell_iterator;
 * const auto cell_worker = [&](const Iterator &cell, ScratchData &scratch,
 *                               CopyData &copy) {
 *   assemble_cell(cell, scratch, copy);
 * };
 * const auto copier = [&](const CopyData &copy) {
 *   constraints.distribute_local_to_global(copy.cell_matrix,
 *                                          copy.local_dof_indices, system_matrix);
 * };
 * WorkStream::run(dof_handler.begin_active(), dof_handler.end(),
 *                 cell_worker, copier, scratch_data, copy_data);
 * @endcode
 *
 * <h3>MeshWorker</h3>
 *
 * MeshWorker::mesh_loop() builds on WorkStream and provides separate callbacks
 * for cells, boundary faces, and interior faces. Flags select which
 * contributions to compute and how interior faces are visited. This is useful
 * for formulations with face terms, such as discontinuous Galerkin methods:
 * the application supplies the local integrals and the loop handles traversal.
 * Scratch and copy data still separate local computation from insertion into
 * global objects. See step-12 and the @ref MeshWorker module.
 *
 * The traversal can include all three kinds of contributions in one call:
 * @code
 * MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
 *                       cell_worker, copy_cell_and_face_data,
 *                       scratch_data, copy_data,
 *                       MeshWorker::assemble_own_cells |
 *                         MeshWorker::assemble_boundary_faces |
 *                         MeshWorker::assemble_own_interior_faces_once,
 *                       boundary_worker, face_worker);
 * @endcode
 * Here the application has also defined @p boundary_worker and @p face_worker
 * with the signatures documented by MeshWorker::mesh_loop(). They append face
 * contributions to @p CopyData; @p copy_cell_and_face_data inserts both cell
 * and face contributions into global objects. This requires additional face
 * storage in the scratch and copy data, as illustrated in step-12.
 *
 * <h3>Matrix-free operator application</h3>
 *
 * When an iterative solver only needs the action of an operator on a vector,
 * it is also possible to avoid assembling and storing the global matrix.
 * MatrixFree organizes cell and face loops and the data needed for them;
 * FEEvaluation and FEFaceEvaluation evaluate finite element functions and
 * integrate local contributions. These contributions are accumulated into a
 * destination vector on each operator application. For tensor-product
 * elements, sum factorization can reduce the cost of local evaluation,
 * particularly at higher polynomial degrees. See step-37 and the
 * @ref matrixfree module.
 *
 * This approach changes the operator representation, not just the assembly
 * loop. It saves matrix storage but requires an operator implementation and
 * a compatible preconditioner. An assembled matrix remains useful when a
 * solver or preconditioner requires its entries. Applications can also
 * combine matrix-free operator application with assembled auxiliary matrices.
 * For the same constant-coefficient Laplace operator, the following member
 * functions of an application class @p LaplaceOperator replace matrix
 * assembly and multiplication. The class has an initialized MatrixFree<dim>
 * member @p data, a compile-time polynomial degree @p fe_degree, and the alias
 * @p VectorType for LinearAlgebra::distributed::Vector<double>. Its vectors
 * use the layout supplied by MatrixFree::initialize_dof_vector().
 * @code
 * void local_apply(const MatrixFree<dim> &mf, VectorType &dst,
 *                  const VectorType &src,
 *                  const std::pair<unsigned int, unsigned int> &range) const
 * {
 *   FEEvaluation<dim, fe_degree> phi(mf);
 *   for (unsigned int cell = range.first; cell < range.second; ++cell)
 *     {
 *       phi.reinit(cell);
 *       phi.gather_evaluate(src, EvaluationFlags::gradients);
 *       for (const unsigned int q : phi.quadrature_point_indices())
 *         phi.submit_gradient(phi.get_gradient(q), q);
 *       phi.integrate_scatter(EvaluationFlags::gradients, dst);
 *     }
 * }
 *
 * void vmult(VectorType &dst, const VectorType &src) const
 * {
 *   data.cell_loop(&LaplaceOperator::local_apply, this, dst, src, true);
 * }
 * @endcode
 * The loop runs over cell batches. The final @p true requests zeroing @p dst
 * before accumulating contributions; cell_loop() also handles the required
 * communication for these distributed vectors. Boundary conditions and
 * constraint treatment must match the operator used by the solver, as in
 * step-37.
 */
