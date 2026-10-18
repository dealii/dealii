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
 * The local finite element calculations can
 * remain the same as in an explicit cell loop. See step-9 and the
 * @ref threads "shared-memory parallelism documentation" for examples and
 * details of the WorkStream model.
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
 */
