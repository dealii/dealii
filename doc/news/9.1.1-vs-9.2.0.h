// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

/**
@page changes_between_9_1_1_and_9_2_0 Changes between Version 9.1.1 and 9.2.0

<p>
This is the list of changes made between the release of deal.II version
9.1.1 and that of 9.2.0. All entries are signed with the names of the
author.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="incompatible"></a>
<h3 style="color:red">Incompatibilities</h3>

<p style="color:red">
Following are a few modifications to the library that unfortunately
are incompatible with previous versions of the library, but which we
deem necessary for the future maintainability of the
library. Unfortunately, some of these changes will require
modifications to application programs. We apologize for the
inconvenience this causes.
</p>

<ol>

 <li>
  Changed: The operator `std::ostream << Vector` was implemented by
  calling Vector::print() with default arguments. This implies a
  particular formatting of the output that was independent of the flags
  set on the output stream (in particular relating to how floating point
  numbers are to be formatted), and also added a newline to the end of
  the output. Both the formatting and the newline are decisions that are
  better made by the user at the location where the output is
  requested. As a consequence, the behavior has now been changed:
  `operator<<` now simply outputs all elements of the vector as
  individual numbers for which the output stream flags are respected,
  and separated by a single space; no endline is appended to the output
  any more.
  <br>
  (Wolfgang Bangerth, 2020/05/05)
 </li>

 <li>
  Removed: The function `MeshWorker::CopyData::operator=(const double)` has been
  removed, as resetting the state of the `copier` from within a `boundary_worker`
  or `face_worker` is not permitted and can lead to hard to trace bugs within a
  user's code.
  <br>
  (Jean-Paul Pelteret, 2020/04/30)
 </li>

 <li>
  Removed: The deprecated functions numbers::is_nan() and
  SLEPcWrappers::SolverBase::set_initial_vector() have been removed.
  <br>
  (Daniel Arndt, 2020/04/24)
 </li>

 <li>
  Removed: The deprecated member functions domain_partitioner(),
  vector_partitioner() and range_partitioner() in TrilinosWrappers classes have
  been removed.
  <br>
  (Daniel Arndt, 2020/04/24)
 </li>

 <li>
  Removed: The deprecated TrilinosWrappers matrix and sparsity patterns and
  reinitialization functions have been removed.
  <br>
  (Daniel Arndt, 2020/04/22)
 </li>

 <li>
  Removed: The deprecated member functions PETScWrappers::VectorBase::ratio() and
  Vector::ratio() have been removed.
  <br>
  (Daniel Arndt, 2020/04/17)
 </li>

 <li>
  Changed: The header file `deal.II/fe/fe_q.h` is not included implicitly via
  `deal.II/fe/mapping_q_generic.h` any more. Add the header
  `deal.II/fe/fe_q.h` if necessary.
  <br>
  (Martin Kronbichler, 2020/04/16)
 </li>

 <li>
  Removed: The deprecated Vector::print overloads have been removed.
  <br>
  (Daniel Arndt, 2020/04/15)
 </li>

 <li>
  Removed: For LinearAlgebra::distributed::Vector, deprecated overloads of the
  member operator=(), sadd(), equ(), local_range(), n_ghost_entries(),
  ghost_elements() and ghost_entries() have been removed.
  <br>
  (Daniel Arndt, 2020/04/13)
 </li>

 <li>
  Removed: The deprecated member functions
  FiniteElement::has_generalized_face_support_points()
  and FiniteElement::get_generalized_face_support_points() have been removed.
  <br>
  (Daniel Arndt, 2020/04/13)
 </li>

 <li>
  Removed: The deprecated functions Utilities::fixed_int_power(),
  LinearAlgebra::distributed::BlockVector::equ(),
  LinearAlgebra::distributed::BlockVector::sadd()
  and FEEvaluationBase::get_cell_data_number() and the deprecated variable
  DoFHandler::invalid_dof_index and hp::DoFHandler::invalid_dof_index have been
  removed.
  <br>
  (Daniel Arndt, 2020/04/08)
 </li>

 <li>
  Changed: The FETools::lexicographic_to_hierarchic_numbering() and
  FETools::hierarchic_to_lexicographic_numbering() functions now take a single
  integer for the degree of the associated FE_Q element and return the
  numbering. The old functions taking the result array as the last argument and
  the ones taking a FiniteElementData as an argument have been deprecated in
  favor of the simpler setup with the degree.
  <br>
  (Martin Kronbichler, 2020/04/07)
 </li>

 <li>
  Removed: The deprecated class PArpackSolver::Shift has been removed.
  Use LinearOperator instead.
  <br>
  (Daniel Arndt, 2020/04/04)
 </li>

 <li>
  Removed: The deprecated classes PointerMatrixBase, PointerMatrix,
  PointerMatrixVector and PointerMatrixAux and the related function
  new_pointer_matrix_base() have been removed. Use LinearOperator instead.
  <br>
  (Daniel Arndt, 2020/04/02)
 </li>

 <li>
  Removed: The deprecated FilteredMatrix class has been removed.
  Use the LinearOperator class instead.
  <br>
  (Daniel Arndt, 2020/04/01)
 </li>

 <li>
  Removed: The deprecated class SwappableVector has been removed.
  <br>
  (Daniel Arndt, 2020/03/31)
 </li>

 <li>
  Removed: The deprecated member function Multigrid::set_debug() has been removed.
  Use the various Multigrid::connect_* functions instead.
  <br>
  (Daniel Arndt, 2020/03/30)
 </li>

 <li>
  Removed: The deprecated classes BlockMatrixArray (and derived classes) and MeanValueFilter
  have been removed. Use BlockLinearOperator instead.
  <br>
  (Daniel Arndt, 2020/03/28)
 </li>

 <li>
  Removed: The deprecated member function
  LAPACKFullMatrix::apply_lu_factorization() has been removed. Use
  LAPACKFullMatrix::solve() instead.
  <br>
  (Daniel Arndt, 2020/03/27)
 </li>

 <li>
  Removed: The deprecated member function Triangulation::set_manifold() taking
  one parameter has been removed. Use Triangulation::reset_manifold() instead.
  <br>
  (Daniel Arndt, 2020/03/27)
 </li>

 <li>
  Removed: The deprecated functions DoFHandler::distribute_mg_dofs() and
  MGTransferPrebuilt::initialize_constraints() and deprecated constructors for
  MGSmootherBlock, MGTransferPrebuilt, MGTransferBlockBase and
  MGTransferBlockSelect have been removed.
  <br>
  (Daniel Arndt, 2020/03/27)
 </li>

 <li>
  Removed: The deprecated member functions Timer::get_data(),
  Timer::get_total_data(), Timer::print_data(), Timer::print_total_data() and
  Timer::operator(), Timer::get_lap_time() have been removed. Use
  Timer::get_last_lap_wall_time_data(), Timer::get_accumulated_wall_time_data(),
  Timer::print_last_lap_wall_time_data(),
  Timer::print_accumulated_wall_time_data() and
  Timer::cpu_time() and Timer::last_wall_time() instead.
  <br>
  (Daniel Arndt, 2020/03/27)
 </li>

 <li>
  Removed: The deprecated method hp::DoFHandler::get_fe() has been removed. Please use
  hp::DoFHandler::get_fe_collection() instead.
  <br>
  (Peter Munch, 2020/03/26)
 </li>

 <li>
  Changed: The function ParameterHandler::parse_input() can now also parse xml- and json-files,
  depending on the file suffix. From now on, only files ending with .prm are treated as standard
  parameter files.
  <br>
  (Niklas Fehn, Peter Munch, 2020/03/26)
 </li>

 <li>
  Removed: The deprecated function Utilities::MPI_Partitioner::get_communicator()
  has been removed. Use Utilities::MPI::Partitioner::get_mpi_communicator()
  instead.
  <br>
  (Daniel Arndt, 2020/03/26)
 </li>

 <li>
  Removed: ParameterHandler::print_parameters_section been removed.
  <br>
  (Daniel Arndt, 2020/03/26)
 </li>

 <li>
  Removed: The deprecated function MGLevelObject::clear() has been removed.
  Use MGLevelObject::clear_elements() instead.
  <br>
  (Daniel Arndt, 2020/03/26)
 </li>

 <li>
  Removed: The deprecated FEValuesBase::get_all_normal_vectors() function has been
  removed. Use FEValuesBase::get_normal_vectors() instead.
  <br>
  (Daniel Arndt, 2020/03/26)
 </li>

 <li>
  Deprecated: SparsityTools::distribute_sparsity_pattern() and
  SparsityTools::gather_sparsity_pattern() taking the locally owned index sets
  of all processes as second arguments have been deprecated. Use the respective
  functions only providing the locally owned dofs of the calling process.
  <br>
  (Martin Kronbichler, 2020/03/22)
 </li>

 <li>
  Changed: VectorizedArray::n_array_elements has been deprecated.
  Please use the method VectorizedArray::size() to access the same information.
  <br>
  (Peter Munch, 2020/03/20)
 </li>

 <li>
  New: The underlying type for types::global_dof_index used to be
  `unsigned long long int`. On all reasonable systems, this is actually
  a 64-bit unsigned integer type, but the name wouldn't tell you that
  unless you already knew. This is now fixed: Instead of the mouthful
  of a type, we now use the more concise `uint64_t`: An unsigned integer
  type with exactly 64 bits of accuracy.
  <br>
  Strictly speaking, this may be an incompatible change: On some
  systems, `unsigned long long int` and `uint64_t` may have the same
  size, but are not the same type. In practice, as long as you have used
  the type `types::global_dof_index`, you will not see any difference at
  all.
  <br>
  (Wolfgang Bangerth, 2020/03/01)
 </li>

 <li>
  Removed: The DataOutBase namespace contained a number of `write_*`
  functions that had already been deprecated in deal.II 9.1. These were
  removed, but all were of mostly internal interest and unlikely to be
  used in user programs.
  <br>
  (Wolfgang Bangerth, 2020/02/05)
 </li>

 <li>
  Changed: FESeries::Fourier and FESeries::Legendre now require to use an
  individual number of coefficients per direction for each finite element
  registered in the provided FECollection.
  <br>
  (Marc Fehling, 2020/01/22)
 </li>

 <li>
  Changed: Adjusted ranges for GridRefinement::refine_and_coarsen_fixed_fraction():
  Now criteria of every cell will be considered.
  <br>
  (Marc Fehling, 12/05/2019)
 </li>

 <li>
  Changed: Adjusted threshold for fixed number adaptation: The last qualifying
  cells for refinement and coarsening now specify their respective threshold.
  This affects both GridRefinement::fixed_number_refinement() and
  parallel::distributed::GridRefinement::fixed_number_refinement().
  <br>
  (Marc Fehling, 12/05/2019)
 </li>

 <li>
  Removed: The type alias FunctionParser::ConstMapIterator has been removed.
  <br>
  (David Wells, 2019/12/05)
 </li>

 <li>
  Deprecated: The include files in `deal.II/base/std_cxx11` have been
  deprecated. Just use the corresponding C++11 header files.
  <br>
  (Wolfgang Bangerth, 2019/08/21)
 </li>

 <li>
  Removed: The include files in `include/deal.II/base/std_cxx1x`,
  deprecated since 2014, have been removed. They had introduced a
  namespace `std_cxx1x` at a time when deal.II did not yet require
  C++11-compatible compilers, and this namespace was later also made
  available via the name `std_cxx11` once C++11 was actually
  published. In any case, none of this is necessary any more since we
  now have C++11 compilers.
  <br>
  (Wolfgang Bangerth, 2019/12/02)
 </li>

 <li>
  Removed: The deprecated method MatrixFree::get_size_info() has been removed.
  Please use the method MatrixFree::get_task_info() to access the same information.
  <br>
  (Peter Munch, 2019/12/01)
 </li>

 <li>
  Removed: The typedefs DataOut::active_cell_iterator and
  DataOut_DoFData::active_cell_iterator have been removed. These typedefs are
  an implementational detail that should have never exposed in the interface.
  <br>
  (Matthias Maier, 2019/11/25)
 </li>

 <li>
  Changed: The variable level_mg_handler has been renamed
  to mg_level in MatrixFree and its associated AdditionalData.
  <br>
  (Peter Munch, 2019/10/17)
 </li>

 <li>
  Changed: The type of exact_solution in VectorTools::integrate_difference() was
  double and fe_function could have any numerical type (double, float,
  std::complex<double> or std::complex<float>). This would lead in certain cases
  to unnecessary casts, for
  example from double to float. In addition it was not
  possible to compare a std::complex exact_solution to a std::complex fe_function.
  Now the types of exact_solution and fe_function must be the same, as a result
  integrate_difference() can be used to compare a
  std::complex exact_solution to
  a std::complex fe_function. The old version of the function has been deprecated.
  <br>
  (Daniel Garcia-Sanchez, 2019/10/09)
 </li>

 <li>
  Changed: Several variations of the
  DoFTools::make_periodicity_constraints() function expressed the
  `direction` as an `int` variable. But it's really an unsigned
  quantity, and as a consequence was changed to `unsigned int`.
  <br>
  (Wolfgang Bangerth, 2019/09/26)
 </li>

 <li>
  Changed: The parallel::TriangulationBase class does not store the number of active cells
  of all MPI processes any more. Instead, the information is computed on-demand when calling
  the function parallel::TriangulationBase::compute_n_locally_owned_active_cells_per_processor.
  <br>
  (Peter Munch, 2019/09/20)
 </li>

 <li>
  Removed: The enum TimestepControl::Strategy is removed.
  <br>
  (Reza Rastak, 2019/09/18)
 </li>

 <li>
  Changed: The
  Algorithms::TimeStepControl::file_name_format() functions have been
  removed. They were setters and getters for an internal field of that
  class that turned out to be unused anywhere -- in other words, the
  functions had no real functionality.
  <br>
  (Wolfgang Bangerth, 2019/09/18)
 </li>

 <li>
  Removed: BlockVectorBase::equ(a, U, b, V) is removed.
  <br>
  (Reza Rastak, 2019/09/03)
 </li>

 <li>
  Changed: The interface to the parallel::CellWeights class was generalized
  by introducing static member functions. This way, users are capable of
  connecting and disconnecting callback functions manually if desired without
  instantiating an object of this class. The legacy interface has been deprecated.
  <br>
  (Marc Fehling, 2019/09/03)
 </li>

 <li>
  Changed: A number of classes implementing tensor-valued polynomials
  had a member function called `compute_n_pols()`. This name has now
  been changed to `n_polynomials()`, better reflecting our usual naming
  scheme as well as what the function does -- which in many cases is not
  actually computing but just returning something.
  <br>
  (Wolfgang Bangerth, 2019/08/2)
 </li>

 <li>
  Changed: The complex operator overloads for multiplication and division are
  restricted to floating point types.
  <br>
  (Daniel Arndt, 2019/08/21)
 </li>

 <li>
  Changed: Various classes implementing tensor polynomials had a
  function `compute()`. This function has been renamed to `evaluate()`
  since this name more adequately represents the fact that the function
  does not in fact compute the polynomial, but only evaluates it for a
  given value of the unknown variable.
  <br>
  (Wolfgang Bangerth, 2019/08/21)
 </li>

 <li>
  Changed: The parallel::Triangulation class has been renamed to parallel::TriangulationBase
  to better reflect its purpose of a base to a series of parallel triangulations such as
  parallel::distributed::Triangulation or parallel::shared::Triangulation.
  <br>
  (Peter Munch, 2019/08/19)
 </li>

 <li>
  Replaced: boost::optional was replaced by std_cxx17::optional that wraps the
  former and std::optional depending on compiler support for C++17.
  <br>
  (Daniel Arndt, 2019/08/17)
 </li>

 <li>
  Changed: The test suite now requires numdiff (i.e., it cannot be run with
  simply diff).
  <br>
  (David Wells, 2019/08/10)
 </li>

 <li>
  Improved: CellDataStorage stores a reference to a single unique triangulation
  to ensure that cells from other triangulations cannot store any data on the
  current object.
  <br>
  (Reza Rastak, 2019/08/08)
 </li>

 <li>
  Changed: mapping.h referred to different kinds of mappings in an enum
  MappingType; however, since it's an enum, MappingType is a misleading name.
  All instances of MappingType and mapping_type have been changed to
  MappingKind and mapping_kind, respectively. The documentation has been
  updated in order to reflect these changes as well.
  <br>
  (Graham Harper, 2019/08/06)
 </li>

 <li>
  Changed: The GridOut::write_svg() function now no longer outputs the
  cell level and number by default; this also obviates the need for the
  colorbar and the caption. This is achieved by changing the defaults in
  the constructor of the GridOutFlags::Svg structure. This likely comes
  closer to what most users want to see from this function -- namely,
  just the mesh for inclusion into text.
  <br>
  (Wolfgang Bangerth, 2019/08/05)
 </li>

 <li>
  Removed: The deprecated PreconditionChebyshev member variables nonzero_starting
  and matrix_diagonal_inverse were removed.
  <br>
  (Daniel Arndt, 2019/07/11)
 </li>

 <li>
  Changed: The CoarseningStrategies struct has been moved out of the
  parallel::distributed::CellDataTransfer class into a separate header and
  is now treated as a namespace. Its static member functions are now free
  functions and only take a `std::vector` as a parameter that contains all
  the data from the children. Therefore, the `coarsening_strategy`
  parameter for the constructor of the
  parallel::distributed::CellDataTransfer class has to be adjusted
  accordingly as well.
  <br>
  (Marc Fehling, 2019/06/26)
 </li>

 <li>
  Removed: The deprecated headers `deal.II/grid/tria_boundary.h` and
  `deal.II/grid/tria_boundary_lib.h` have been removed.
  <br>
  (Martin Kronbichler, 2019/06/13)
 </li>

 <li>
  Fixed: During coarsening, we need to decide which finite element is assigned
  to parent cells. The decision criterion was not consistent, but is now unified:
  The least dominant element of their children will be picked, according to
  hp::FECollection::find_dominated_fe_extended(). This affects the implementations
  of SolutionTransfer, parallel::distributed::SolutionTransfer, parallel::CellWeights,
  and hp::DoFHandler.
  <br>
  (Marc Fehling, 2019/06/06)
 </li>

 <li>
  Changed: The functions DoFHandler::n_locally_owned_dofs_per_processor(),
  DoFHandler::locally_owned_dofs_per_processor() and
  DoFHandler::locally_owned_mg_dofs_per_processor() previously returned a
  reference to an internally stored array of index sets on all processors. As
  this cannot scale to large processor counts, these functions have been marked
  deprecated and only populate the internal vectors on the first demand. This
  means that the first call must be a collective call on all MPI ranks to ensure
  that the underlying MPI_Allgather function completes. Use the
  new functions DoFHandler::compute_n_locally_owned_dofs_per_processor(),
  DoFHandler::compute_locally_owned_dofs_per_processor() and
  DoFHandler::compute_locally_owned_mg_dofs_per_processor() instead or, even
  better for scalability, avoid them in favor of some local communication.
  <br>
  (Martin Kronbichler, 2019/06/03)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>

 <li>
  New: The new Differentiation::SD::BatchOptimizer class can be used to accelerate
  (in some cases, significantly) evaluation of the symbolic expressions using an
  assortment of techniques. At the moment, this includes using common
  subexpression elimination to prevent, as much as possible, repetitive evaluation
  of subexpressions that are found in one or more expressions. The expressions
  may also be transformed into other equivalent data structures that are simply
  less computationally costly than their original symbolic expression tree.
  It is also possible to compile a set of expressions using the LLVM JIT compiler,
  rendering near-native evaluation performance at the cost of the compilation
  itself.
  <br>
  (Jean-Paul Pelteret, 2020/05/02)
 </li>

 <li>
  New: The Step-50 tutorial has been updated. We discuss how to
  use the multilevel preconditioner from step-16 in parallel and
  give a comparison of parallel scaling between GMG and AMG for
  a 3D Laplace example where adaptive refinement comes from a
  residual-based, cell-wise a posteriori error estimator.
  <br>
  (Thomas C. Clevenger, Timo Heister, 2020/04/20)
 </li>

 <li>
  New: The step-58 tutorial program demonstrates the solution of the
  nonlinear Schr&ouml;dinger equation.
  <br>
  (Wolfgang Bangerth, Yong-Yong Cai, 2020/04/15)
 </li>

 <li>
  New: The DataOut class has learned to deal with complex-valued
  solution vectors in all regards. This includes writing the real and
  imaginary parts of not only scalar complex-valued fields, but also
  complex-valued vector and tensor fields. DataOut can now also pass
  such fields to classes derived from DataPostprocessor.
  <br>
  (Wolfgang Bangerth, 2020/04/03)
 </li>

 <li>
  New: The step-53 tutorial in the form of a Jupyter notebook
  using deal.II Python interface.
  <br>
  (Alexander Grayver, 2020/03/30)
 </li>

 <li>
  New: The step-67 tutorial program presents an explicit time integrator for the
  compressible Euler equations discretized with a high-order discontinuous
  Galerkin scheme using the matrix-free infrastructure. Besides the use of
  matrix-free evaluators for systems of equations and over-integration, it also
  presents MatrixFreeOperators::CellwiseInverseMassMatrix, a fast implementation
  of the action of the inverse mass matrix in the DG setting using tensor
  products.
  <br>
  (Martin Kronbichler, 2020/03/28)
 </li>

 <li>
  Removed: The outdated \step-19 that was used to combine parallel visualization
  output after a simulation has been removed. The DataOut class can be used for
  parallel output.
  <br>
  (Timo Heister, 2020/03/18)
 </li>

 <li>
  New: The step-69 tutorial program presents a first-order accurate
  <i>guaranteed maximum wavespeed method</i> based on a first-order <i>graph
  viscosity</i> for solving Euler's equations of gas dynamics. It introduces
  and discusses number of interesting programming techniques: hybrid
  parallelization (MPI + task based thread parallelization); the use of
  non-distributed deal.II matrices in an MPI parallel context with index
  rewriting; offline data preprocessing; background-thread offloading for
  output; and checkpointing/restore.
  <br>
  (Matthias Maier, Ignacio Tomas, 2020/02/19)
 </li>

 <li>
  Fixed: Large computations with more than 4 Billion cells or DoFs
  would run into issues in various situations. Bugs in BlockVector,
  Triangulation, Multigrid, GridRefinement, and other places are
  now fixed.
  <br>
  (Thomas Clevenger, Timo Heister, Martin Kronbichler, Peter Munch, 2020/02/19)
 </li>

 <li>
  New: The step-49 tutorial in the form of a Jupyter notebook
  using deal.II Python interface.
  <br>
  (Alexander Grayver, 2020/02/12)
 </li>

 <li>
  New: The step-47 tutorial program demonstrates the solution of the
  biharmonic equation, a fourth-order differential equation.
  <br>
  (Natasha Sharma, Guido Kanschat, Timo Heister, Wolfgang Bangerth, Zhuoran Wang, 2020/01/15)
 </li>

 <li>
  New: The SymbolicFunction<dim> class allows one to leverage the SymEngine library to generate
  dealii::Function objects where the gradients, Laplacians, and Hessians are computed symbolically
  providing also the possibility to extract the time derivative of the SymbolicFunction<dim> object
  as another SymbolicFunction<dim> object.
  <br>
  (Luca Heltai, 2019/12/14)
 </li>

 <li>
  New: Implementation of the IDR(s) Krylov solver for non-symmetric,
  indefinite linear systems.
  <br>
  (Thomas Clevenger, 2019/12/05)
 </li>

 <li>
  New: Particles::ParticleHandler::insert_global_particles() allows one to pass to a ParticleHandler object a vector
  of positions and a vector of properties to insert. Differently from the other `insert*` methods,
  this one allows particles to fall within artificial cells. In this case, the method infers who
  should receive the positions and the properties, and sends this information to that process.
  This is useful when constructing particles from non-matching triangulations, where the
  distribution of the particle positions is arbitrary, and possibly not related to the locally owned cells.
  <br>
  (Luca Heltai, Bruno Blais, 2019/11/25)
 </li>

 <li>
  New: Add new GridGenerator tool to create a C-type mesh around a NACA or
  Joukowski airfoil.
  <br>
  (Elias Dejene, Peter Munch, 2019/11/12)
 </li>

 <li>
  New: Namespace SmoothnessEstimator providing smoothness estimation strategies
  for hp-adaptive FEM based on Fourier and Legendre series expansions.
  <br>
  (Marc Fehling, Denis Davydov, 2019/10/30)
 </li>

 <li>
  New: The step-12 tutorial program has been changed to use
  FEInterfaceValues. The old version of step-12 is still available as step-12b.
  <br>
  (Timo Heister, 2019/08/27)
 </li>

 <li>
  New: The FEInterfaceValues class provides a new abstraction to assemble
  interface terms between two neighboring cells. This is commonly used in
  Discontinous Galerkin methods.
  <br>
  (Timo Heister, 2019/08/24)
 </li>

 <li>
  New: Introduce a new distributed triangulation class parallel::fullydistributed::Triangulation.
  It partitions the coarse-grid as well as works for hanging nodes, geometric multigrid, and
  arbitrary dimensions.
  <br>
  (Peter Munch, 2019/07/26)
 </li>

 <li>
  Improved: deal.II now bundles a subset of Boost 1.70 instead of a subset
  of Boost 1.62.
  <br>
  (Daniel Arndt, 2019/07/19)
 </li>

 <li>
  New: The step-65 tutorial program presents TransfiniteInterpolationManifold, a
  manifold class that can propagate curved boundary information into the
  interior of a computational domain, and MappingQCache for fast operations for
  expensive manifolds.
  <br>
  (Martin Kronbichler, 2019/06/06)
 </li>

 <li>
  New: The namespace Utilities::MPI::ConsensusAlgorithms has been added. It
  provides efficient implementations for communication patterns (PEX and NBX) to
  retrieve data from other processes in a dynamic-sparse way.
  <br>
  (Peter Munch, Martin Kronbichler, 2019/06/03)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>

 <li>
  Fixed: GridTools::find_closest_vertex_of_cell() now takes an optional mapping argument, to take into
  account of all those cases where a mapping may modify the location of the cell vertices.
  <br>
  (Luca Heltai, 2020/05/10)
 </li>

 <li>
  Fixed: Fix a bug in tutorial step-35, which caused incorrect initial condition initialization. This bug could lead to incorrect results if the user decided to edit the equation data.
  <br>
  (Alexey Ozeritskiy, 2020/05/09)
 </li>

 <li>
  Fixed: FEFaceEvaluation::read_dof_values_plain() would sometimes access
  invalid dof indices and crash or produce wrong results for cells with
  constraints applied. This is now fixed.
  <br> (Martin Kronbichler, Niklas Fehn, 2020/05/07)
 </li>

 <li>
  Changed: Allow hp::DoFHandler to distribute dofs based on a triangulation class derived from parallel::DistributedTriangulationBase other than parallel::distributed::Triangulation.
  <br>
  (Sebastian Stark, 2020/05/06)
 </li>

 <li>
  New: Add a flag to overlap MPI communication and computation when using CUDA
  matrix-free and CUDA-aware MPI.
  <br>
  (Peter Munch and Bruno Turcksin, 2020/05/06)
 </li>

 <li>
  Changed: The order by which MatrixFree::cell_loop() and MatrixFree::loop()
  pass through cells has been changed to preferably group cells with the same
  parent into the same batch of cells with vectorization, which increases data
  locality and slightly improves performance.
  <br>
  (Martin Kronbichler, 2020/05/06)
 </li>

 <li>
  New: Codim one instantiations for make_hanging_node_constraints.
  <br>
  (Sebastian Stark, 2020/05/05)
 </li>

 <li>
  Changed: Compute middle vertices during mesh refinement of Triangulation<2,3>
  using transfinite interpolation to be consistent with what happens to quads
  for Triangulation<3,3>
  <br>
  (Sebastian Stark, 2020/05/04)
 </li>

 <li>
  Improved: MatrixFreeOperators::LaplaceOperator now also supports cell-wise
  constant coefficient, which can be used by a Table with a single column
  per cell.
  <br>
  (Martin Kronbichler, 2020/05/03)
 </li>

 <li>
  New: ExtractLevelVisitor and extract_rtree_level() allow one to return the vector of BoundingBox objects
  that make up a specific level of a boost::rtree object.
  <br>
  (Luca Heltai, 2020/05/03)
 </li>

 <li>
  New: There is a new triad of functions FEValuesBase::dof_indices(),
  FEValuesBase::dof_indices_starting_at() and
  FEValuesBase::dof_indices_ending_at()
  that makes writing range-based for loops over all local cell degrees of
  freedom much simpler.
  <br>
  (Jean-Paul Pelteret, 2020/04/29)
 </li>

 <li>
  New: There is a new function FEValuesBase::quadrature_point_indices()
  that makes writing range-based for loops over all quadrature points of a
  cell or face much simpler.
  <br>
  (Jean-Paul Pelteret, 2020/04/29)
 </li>

 <li>
  New: Add Utilities::create_evenly_distributed_partitioning(my_partition_id, n_partitions, total_size)
  to create a one-to-one evenly distributed ascending partitioning from total size.
  Add Utilities::MPI::create_evenly_distributed_partitioning(comm, total_size) to use processor ID and
  number of MPI processes to determine partitioning.
  <br>
  (Doug Shi-Dong, 2020/04/29)
 </li>

 <li>
  Fixed: Copying objects of type hp::FEValuesBase and derived classes
  led to unexpected results because both source and destination kept
  pointers to a shared space. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/04/28)
 </li>

 <li>
  New: All FEValues objects in a hp::FEValues collection can now be
  pre-computed with hp::FEValues::precalculate_fe_values().
  <br>
  (Marc Fehling, 2020/04/27)
 </li>

 <li>
  New: ParameterAcceptor::initialize() now uses ParameterHandler::print_parameters() internally.
  <br>
  (Niklas Fehn, Luca Heltai, 2020/04/24)
 </li>

 <li>
  Improved: Add Tvmult operation to PETScWrappers::PreconditionerBase.
  <br>
  (Daniel Garcia-Sanchez, 2020/04/22)
 </li>

 <li>
  Improved: GridGenerator::hyper_shell() in 3d now supports more `n_cells`
  options. While previously only 6, 12, or 96 cells were possible, the function
  now supports any number of the kind $6 \times 2^m$ with $m$ a non-negative
  integer. The new cases $m=2,3$ and $m\geq 5$ correspond to refinement in the
  azimuthal direction of the 6 or 12 cell case with a single mesh layer in
  radial direction, and are intended for shells that are thin and should be
  given more resolution in azimuthal direction.
  <br>
  (Martin Kronbichler, 2020/04/07)
 </li>

 <li>
  Fixed: VectorTools::get_position_vector would return wrong positions for
  curved cells in the interior of the domain. This is now fixed.
  <br> (Martin Kronbichler, Doug Shi-Dong, 2020/04/05)
 </li>

 <li>
  New: Add new function Utilities::MPI::compute_set_union(), which computes the union of given input
  sets and vectors of all processes in a given MPI communicator.
  <br>
  (Peter Munch, 2020/04/04)
 </li>

 <li>
  New: There is a new function GridGenerator::hyper_ball_balanced() which uses a
  more refined mesh of 12 cells in 2D and 32 in 3D to create a disk in 2D and
  ball in 3D, producing better mesh quality than the default
  GridGenerator::hyper_ball().
  <br> (Martin Kronbichler, 2020/04/02)
 </li>

 <li>
  Changed: The vertex positions of GridGenerator::quarter_hyper_ball() have been
  slightly changed. It concerns the points along the straight faces and in the
  interior. The new positions have been determined to balance the sizes and
  aspect ratios (expressed by the minimal and maximal singular values of the
  Jacobian from reference to real coordinates), producing a slightly better
  mesh quality.
  <br> (Martin Kronbichler, 2020/04/02)
 </li>

 <li>
  New: Utilities::encode_base64() and Utilities::decode_base64() allow to encode and decode binary strings
  to and from Base64 ASCII format.
  <br>
  (Luca Heltai, 2020/04/01)
 </li>

 <li>
  New: Added Particles::Utilities::create_interpolation_sparsity_pattern() and
  Particles::Utilities::create_interpolation_matrix() that create interpolation sparsities and
  matrices between ParticleHandler and DoFHandler objects.
  <br>
  (Bruno Blais, Luca Heltai, 2020/03/31)
 </li>

 <li>
  Fixed: A bug in the setup of CellAccessor::level_subdomain_id() for multigrid
  levels that manifested itself only on non-refined meshes has been fixed.
  <br> (Martin Kronbichler, 2020/03/30)
 </li>

 <li>
  Fixed: FEValues would previously detect cell similarities when combined with
  MappingFEField and MappingManifold when the underlying mesh had similarities,
  even though the deformed cells did not. This would lead to wrong Jacobians and
  transformed shape function gradients, among others, when running deal.II
  without threads. This is now fixed.
  <br> (Martin Kronbichler, 2020/03/30)
 </li>

 <li>
  Fixed: step-26 and step-52 now also output the simulation time
  in the VTK/VTU output files.
  <br>
  (Wolfgang Bangerth, 2020/03/27)
 </li>

 <li>
  Improved: Step-21 has been adapted to use the new class DiscreteTime.
  <br>
  (Reza Rastak, 2020/03/26)
 </li>

 <li>
  Extension: An optional additional parameter assert_mandatory_entries_are_found
  is introduced to the function ParameterHandler::parse_input() to avoid errors,
  e.g., due to typos in the input file.
  <br>
  (Niklas Fehn, 2020/03/26)
 </li>

 <li>
  Changed: The SIMD vectorization capabilities in deal.II are now controlled by
  the macro DEAL_II_VECTORIZATION_WIDTH_IN_BITS containing the numbers of bits
  in the SIMD arrays (e.g., 256 for AVX or 512 for AVX-512).
  <br> (Martin Kronbichler, 2020/03/24)
 </li>

 <li>
  New: PreconditionChebychev can now report information about eigenvalue and
  degree estimation when calling PreconditionChebychev::estimate_eigenvalues().
  <br> (Timo Heister, 2020/03/23)
 </li>

 <li>
  New: ParameterHandler: add functionality ensuring that parameters are
  indeed set in order to avoid errors that remain unrecognized otherwise.
  By extending the interfaces of declare_entry() and add_parameter(),
  parameters can be declared as mandatory parameters that must be set
  when parsing from an input file or by a set() function.
  <br>
  (Niklas Fehn, 2020/03/23)
 </li>

 <li>
  New: The function ParameterHandler::print_parameters() can now also
  print a reduced parameter tree for XML and JSON file formats.
  <br>
  (Peter Munch, 2020/03/22)
 </li>

 <li>
  New: Variants of SparsityTools::distribute_sparsity_pattern() and
  SparsityTools::gather_sparsity_pattern() taking only the locally owned index
  set of the calling MPI process have been added. They avoid the previously
  necessary MPI all-to-all communication to compute the owned dofs on all
  processes and instead go through Utilities::MPI::compute_index_owned with
  sparse communication. On more than around 1,000 MPI ranks, the new functions
  should be considerably faster.
  <br>
  (Martin Kronbichler, 2020/03/22)
 </li>

 <li>
  Added high-order mesh output to the step-11 tutorial.
  <br>
  (Alexander Grayver, 2020/03/20)
 </li>

 <li>
  New: DataOut::add_mg_data_vector() allows the user to output solutions on
  multigrid levels.
  <br> (Timo Heister, 2020/03/18)
 </li>

 <li>
  Fixed: Add missing fields to MatrixFree::copy_from.
  <br>
  (Daniel Jodlbauer, 2020/03/16)
 </li>

 <li>
  Improved: The detection of compression in the mapping data structures of
  MatrixFree has been made more general. It can now also detect when two cells
  (not necessarily adjacent ones) are translations of each other, so that the
  same set of Jacobians can be used on both. This helps the performance of
  matrix-free operator evaluation in case the geometry fits into processor
  caches and the operation is memory-bandwidth limited.
  <br>
  (Martin Kronbichler, 2020/03/14)
 </li>

 <li>
  Improved: The initialization of the mapping part of MatrixFree for mappings of
  type MappingQGeneric or derived classes has been made much faster. Apart from
  a single query per cell for the geometry via the MappingQGeneric class, all
  interpolations to metric terms at quadrature points are now done with
  matrix-free routines using sum factorization, providing similar performance as
  an operator evaluation.
  <br>
  (Martin Kronbichler, 2020/03/14)
 </li>

 <li>
  Improved: ParticleHandler to have a GridTools::Cache to ensure better efficiency
  <br>
  (Bruno Blais, 2020/03/12)
 </li>

 <li>
  New: LinearAlgebra::set_zero_mean_value()
  Add free function that allows to shift a vector by a constant value in a
  way that its mean value becomes zero. This function only makes use of the
  pure virtual functions add() and mean_value() of the abstract class
  VectorSpaceVector.
  <br>
  (Niklas Fehn, 2020/03/06)
 </li>

 <li>
  New: Introduce a new communication-pattern class Utilities::MPI::NoncontiguousPartitioner.
  It is similar to Utilities::MPI::Partitioner, however, does not make any
  restrictions regarding to the ordering of the underlying index sets. This class enables
  efficient repartitioning of vectors, which might be beneficial in interfacing with
  external libraries that expect a certain fixed order (like checkerboard partitioning).
  <br>
  (Peter Munch, 2020/03/04)
 </li>

 <li>
  New: There are new functions GeometryInfo::face_indices() and
  GeometryInfo::vertex_indices() that make writing range-based for loops
  over all faces or vertices of a cell much simpler.
  <br>
  (Wolfgang Bangerth, 2020/03/04)
 </li>

 <li>
  Fixed: FEFaceEvaluation::read_dof_values() and
  FEFaceEvaluation::distribute_local_to_global() would not correctly represent
  cases with boundary integrals and constraints on the same boundary. This is
  now fixed.
  <br>
  (Niklas Fehn, Martin Kronbichler, 2020/03/03)
 </li>

 <li>
  New: There is now a new type types::global_cell_index that is used to
  denote the global index of a cell in a parallel triangulation, and
  consequently also places where we compute or return the *total* number
  of cells.
  <br>
  (Wolfgang Bangerth, 2020/03/01)
 </li>

 <li>
  Changed: Fixed VectorTools::project_boundary_values_curl_conforming_l2 to
  work with FESystems containing an FE_Nedelec.
  <br>
  (Winnifried Wollner, 2020/02/28)
 </li>

 <li>
  New: GeometryInfo<dim> now holds an array
  GeometryInfo<dim>::unit_tangential_vectors with the tangential vectors on the
  reference cell.
  <br>
  (Martin Kronbichler, 2020/02/28)
 </li>

 <li>
  New: There is now a new function MatrixFree::update_mapping() that allows to
  refresh the stored geometry data when the mapping has changed, for example
  through a Eulerian mesh motion, while keeping all other data structures such
  as vector exchange pattern or dof indices alive.
  <br>
  (Martin Kronbichler, 2020/02/25)
 </li>

 <li>
  New: An IndexSet::split_by_block() function is added which partitions the set
  indices, for example based on a block structure, into blocks.
  <br>
  (Katrin Mang, 2020/02/14)
 </li>

 <li>
  Add: python binds for Triangulation::create_triangulation, GridTools::scale
  and GridGenerator::hyper_cube_with_cylindrical_hole.
  <br>
  (Alexander Grayver, 2020/02/11)
 </li>

 <li>
  Fixed: The output message generated by ExcInvalidBoundaryFunction() in deal.II/numerics/error_estimator.h, i.e. "However, the finite element in use has arg2 components" should be "However, the finite element in use has arg3 components". This is now fixed.
  <br>
  (Jihuan Tian, 2020/02/09)
 </li>

 <li>
  Improved: Add a constructor to ArrayView that takes as argument
  a reference to an arbitrary value type (like int, double, ...).
  <br>
  (Peter Munch, 2020/01/02)
 </li>

 <li>
  New: The function Utilities::MPI::min_max_avg() can now compute the sum, average, minimum, maximum, and the
  process id of the minimum and maximum for each entry of a ArrayView.
  <br>
  (Peter Munch, 2020/02/03)
 </li>

 <li>
  Changed: The GridOutFlags::Svg constructor may now take an additional
  argument supporting the output of boundary ids on faces.
  <br>
  (Graham Harper, 2020/02/02)
 </li>

 <li>
  Fixed: Calling DataOut::merge() with objects that were created using a
  mapping and more than one number of subdivisions greater than one
  sometimes resulted in an unwarranted triggering of an assertion. This
  is now fixed.
  <br>
  (Paras Kumar, Wolfgang Bangerth, 2020/02/02)
 </li>

 <li>
  New: added ParameterHandler::subsection_path_exists() method.
  <br>
  (Pasquale Claudio Africa, 2020/01/29)
 </li>

 <li>
  Add an optional std::function parameter to make_flux_sparsity_pattern which can
  be used to specify over which faces there should be a flux coupling in the
  created sparsity pattern.
  <br>
  (Simon Sticko, 2020/01/28)
 </li>

 <li>
  Improved: ParameterHandler class can now optionally print parameters in the same
  order as they are declared, rather than sorted in alphabetical order.
  <br>
  (Pasquale Claudio Africa, 2020/01/28)
 </li>

 <li>
  Improved: ParameterHandler::parse_input_from_json() can handle values of type `"a" : 1', in addition to
  `"a" : {"value" : 1}'.
  <br>
  (Peter Munch, 2020/01/27)
 </li>

 <li>
  New: ParameterHandler::parse_input_from_json() and  ParameterHandler::parse_input_from_xml() can skip
  undefined sections and entries.
  <br>
  (Peter Munch, 2020/01/27)
 </li>

 <li>
  New: Python bindings for the Quadrature class.
  <br>
  (Alexander Grayver, 2020/01/20)
 </li>

 <li>
  New: Add method to compute measure of distorted quads embedded in a three dimensional space.
  <br>
  (Nicola Giuliani, 2020/01/23)
 </li>

 <li>
  New: FESeries::process_coefficients() now allows to ignore coefficients
  below an absolute threshold, which has to be provided as a parameter.
  <br>
  (Marc Fehling, 2020/01/21)
 </li>

 <li>
  New: GridOut::write_vtu() and GridIn::read_vtu() allow to save and restore a locally refined triangulation,
  using a xml section in the vtu file that is ignored by vtu readers.
  <br>
  (Luca Heltai, Nicola Giuliani, 2020/01/17)
 </li>

 <li>
  New: Utilities::compress() and Utilities::decompress() allow to compress and decompress a string using gzip.
  <br>
  (Luca Heltai, Nicola Giuliani, 2020/01/17)
 </li>

 <li>
  New: FESeries::Fourier and FESeries::Legendre have been extended by an
  equality operator (==) and the functionality to save and load previously
  calculated transformation matrices.
  <br>
  (Marc Fehling, 2020/01/16)
 </li>

 <li>
  New: Equality operator (==) for hp::QCollection.
  <br>
  (Marc Fehling, 2020/01/16)
 </li>

 <li>
  Deprecated: The CellAccessor::active() function did not satisfy the
  usual naming scheme of that class in which functions are called
  CellAccessor::is_ghost(), CellAccessor::is_locally_owned(), etc. As a
  consequence, there is now a new function CellAccessor::is_active(),
  with the old function deprecated.
  <br>
  (Wolfgang Bangerth, 2020/01/16)
 </li>

 <li>
  Fixed: The GridOut::write_svg() function was not instantiated for the
  case `dim=3, spacedim=3` by accident. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/01/14)
 </li>

 <li>
  New: Added a new function GridGenerator::replicate_triangulation() for creating
  a new triangulation by copying a given triangulation repeatedly along the
  coordinate axes.
  <br>
  (David Wells, Victor Zheng, 2020/01/13)
 </li>

 <li>
  Fixed: DataOut::write_vtu_in_parallel() used to generate invalid .vtu files if
  some processors had 0 patches to write. This is now fixed.
  <br> (Timo Heister, 2020/01/12)
 </li>

 <li>
  Improved: GridTools::delete_duplicated_vertices() now runs, for cubelike
  geometries, in $O(n^{3/2})$ time in 2D and $O(n^(5/3))$ time in 3D instead
  of $O(n^2)$ time.
  <br>
  (David Wells, 2020/01/12)
 </li>

 <li>
  New: The mechanism implemented by DataOut via the `first_cell()` and
  `next_cell()` virtual functions has been deprecated. In its place, the
  DataOut::set_cell_selection() function allows doing the same but
  without the need for deriving classes.
  <br>
  (Wolfgang Bangerth, 2020/01/11)
 </li>

 <li>
  New: ParticleHandler::get_particle_positions() and ParticleHandler::set_particle_positions() allow
  to set and get particle positions from various types of sources.
  <br>
  (Bruno Blais, Luca Heltai, 2020/01/10)
 </li>

 <li>
  New: There are new versions of DoFTools::count_dofs_per_component()
  and DoFTools::count_dofs_per_block() that return their results as a
  `std::vector` object, rather than doing so through an output
  argument. The versions that did the latter have been deprecated.
  <br>
  (Wolfgang Bangerth, 2020/01/09)
 </li>

 <li>
  New: A variant of CellwiseInverseMassMatrix::apply() defining the inverse of
  the diagonal by FEEvaluationBase::JxW() from the provided FEEvaluationBase
  variable has been implemented.
  <br>
  (Martin Kronbichler, 2020/01/05)
 </li>

 <li>
  New: Add method is_ancestor_of to CellId.
  <br>
  (Peter Munch, 2020/01/02)
 </li>

 <li>
  New: The function MatrixFree::get_dof_handler() has been templated to enable
  to return either a DoFHandler<dim> or an hp::DoFHandler<dim>.
  <br>
  (Bruno Turcksin, 2020/01/02)
 </li>

 <li>
  New: The function Threads::TaskGroup::size() allows querying the size
  of a task group.
  <br>
  (Wolfgang Bangerth, 2019/12/31)
 </li>

 <li>
  Improved: mu_parser functions now use random number generation facilities
  provided by the standard library.
  <br>
  (Reza Rastak, 2019/20/30)
 </li>

 <li>
  New: The SparseDirectUMFPACK class is now able to solve complex-valued
  problems as well.
  <br>
  (Wolfgang Bangerth, 2019/12/30)
 </li>

 <li>
  New: The method memory_consumption has been implemented in and added
  to a variety of classes, e.g. FE_DGQ.
  <br>
  (Peter Munch, 2019/12/22)
 </li>

 <li>
  New: There are new versions of DoFTools::extract_dofs() that return
  their results as an IndexSet, rather than doing so through an output
  argument of type `std::vector<bool>`. The versions that did the latter
  have been deprecated.
  <br>
  (Wolfgang Bangerth, 2019/12/19)
 </li>

 <li>
  New: The TensorFunctionParser class allows to read
  in a function similarly to the FunctionParser class using the MuParser.
  <br>
  (Konrad Simon, 2019/12/17)
 </li>

 <li>
  Added python wrappers for the FunctionManifold class
  <br>
  (Alexander Grayver, 2019/12/17)
 </li>

 <li>
  Added: Python wrappers for the CellAccessor's methods: active, level,
  index, vertex_index, neighbor_is_coarser, neighbor_of_neighbor and
  GridTools' methods: find_cells_adjacent_to_vertex, minimal_cell_diameter,
  maximal_cell_diameter, and Mapping::project_real_point_to_unit_point_on_face
  <br>
  (Alexander Grayver, 2020/01/06)
 </li>

 <li>
  Improved: GridGenerator::channel_with_cylinder() now places some of the points
  at the transition between the cylinder mesh and the outer structured mesh
  slightly farther away from the cylinder center, thereby improving mesh
  quality.
  <br>
  (Martin Kronbichler, 2019/12/16)
 </li>

 <li>
  New: IndexSet::tensor_product() creates a new IndexSet with global size equal to
  this->size()*other.size(), containing for every element n of this IndexSet, the indices
  of the other IndexSet, contained in the interval [n*other.size(), (n+1)*other.size()).;<br>;(Luca Heltai, 2019/12/12)
 </li>

 <li>
  Fixed: The introduction of step-18 had a mistake in the description of
  the weak formulation. However, given that we choose the boundary forces
  equal to zero, the program correctly implements the equations.
  <br>
  (Ming Yang, 2019/12/06)
 </li>

 <li>
  Added: Python wrappers for the CellAccessor::center, GridTools::transform,
  GridTools::distort_random, GridGenerator::extrude_triangulation,
  GridTools::find_active_cell_around_point
  <br>
  (Alexander Grayver, 2019/12/06)
 </li>

 <li>
  Improved: GridGenerator::subdivided_hyper_cube() gained an option to also set boundary ids.
  <br>
  (Bruno Blais, 2019/12/05)
 </li>

 <li>
  Changed: step-61 no longer generates two output files for the interior
  pressures and the Darcy velocities, respectively. Instead, it puts all
  of this information into a single graphical output file.
  <br>
  (Wolfgang Bangerth, 2019/12/02)
 </li>

 <li>
  Changed: The step-61 program now has a separate function
  `compute_postprocessed_velocity()` that computes the velocity field as
  a postprocessing step after computing the pressure variable. The
  resulting velocity field is used in both computing errors and when
  creating graphical output.
  <br>
  (Wolfgang Bangerth, 2019/12/02)
 </li>

 <li>
  New: Particles::Particle::set_id() allows to set a particle ID number.
  <br>
  (Luca Heltai, 2019/12/02)
 </li>

 <li>
  Fixed: Utilities::MPI::some_to_some() now allows sending to our own mpi process.
  <br>
  (Luca Heltai, 2019/12/02)
 </li>

 <li>
  Fixed: ParticleHandler::insert_particles() taking a std::vector of Point objects
  starts numbering the particles with 0 instead of 1 to be consistent with the
  other ParticleHandler::insert_particles() function and the
  ParticleHandler::insert_particle() function.
  <br>
  (Daniel Arndt, 2019/11/30)
 </li>

 <li>
  Improved: Allow using ParticleHandler with serial triangulations.
  <br>
  (Daniel Arndt, 2019/11/29)
 </li>

 <li>
  Added: New python wrappers for the Manifold class
  with a possibility to assign manifold object to a
  python triangulation object.
  <br>
  (Alexander Grayver, 2019/11/28)
 </li>

 <li>
  Improved: The setup of MappingQGeneric::InternalData within the constructor of
  FEValues would previously unconditionally allocate memory for all shape
  functions and all quadrature points, also for the case where we use the tensor
  product and the full interpolation is unnecessary. This has been fixed,
  improving the situation for very high orders and numbers of quadrature points
  (e.g., avoiding 400 MB of memory for mapping degrees of 15 with $16^3$
  quadrature points in 3D).
  <br>
  (Martin Kronbichler, 2019/11/26)
 </li>

 <li>
  Fixed: MatrixFree::reinit() would sometimes perform invalid index accesses on
  continuous elements with hanging nodes when data structures for face integrals
  are requested. This is now fixed.
  <br>
  (Martin Kronbichler, Peter Munch, Laura Prieto Saavedra, 2019/11/26)
 </li>

 <li>
  New: ParticleHandler::locally_relevant_ids() generates an IndexSet of the locally owned
  particles. This can be used to construct linear algebra vectors to be used with particles.
  <br>
  (Luca Heltai, Bruno Blais, 2019/11/18)
 </li>

 <li>
  Added: New python wrappers for the MappingQGeneric class
  <br>
  (Alexander Grayver, 2019/11/18)
 </li>

 <li>
  Changed: Indent scripts now also format files related to the python bindings
  <br>
  (Alexander Grayver, 2019/11/15)
 </li>

 <li>
  Fixed: MappingQGeneric::transform_real_to_unit_cell() would previously ignore
  the moved vertex locations of derived classes such as MappingQEulerian for
  linear mappings in 1D and 2D. This is now fixed.
  <br>
  (Martin Kronbichler, 2019/11/14)
 </li>

 <li>
  New: DoFTools::locally_owned_dofs_per_component() returns a vector of IndexSet, containing the locally owned
  dofs that refers to the ComponentMask passed as an argument.
  <br>
  (Bruno Blais, Luca Heltai, 2019/11/26)
 </li>

 <li>
  Added: extend python bindings to allow for working with cell's
  neighbors, iterate over faces and modify boundary ids.
  <br>
  (Alexander Grayver, 2019/11/14)
 </li>

 <li>
  New: Added reading writing capabilities for STL files.
  <br>
  (Nicola Giuliani, 2019/11/13)
 </li>

 <li>
  New: add more quadrature choices to QuadratureSelector class.
  <br>
  (Daniel Jodlbauer, 2019/11/11)
 </li>

 <li>
  Fixed: In Step-54, the std namespace was missing for some calls to cout and
  endl.
  <br>
  (Mohammed Hassan, Bruno Turcksin, 2019/11/07)
 </li>

 <li>
  New: template argument for type of vectorized array (vector width) in matrix-free operators.
  <br>
  (Daniel Jodlbauer, 2019/11/06)
 </li>

 <li>
  Fixed: Computations with block vectors with more than 4 billion degrees of
  freedom reported incorrect values for locally_owned_elements(). This is now
  fixed by implementing IndexSet::add_indices() for large offsets.
  <br>
  (Timo Heister, 2019/11/05)
 </li>

 <li>
  New: The Utilities::MPI::CollectiveMutex helper class allows
  protection of critical sections in code using MPI.
  <br>
  (Timo Heister, 2019/11/04)
 </li>

 <li>
  Fixed: FESubfaceValues::reinit() would previously throw an exception when
  called on a face behind a refined periodic neighbor. This case is now
  handled properly.
  <br>
  (Peter Munch, 2019/11/03)
 </li>

 <li>
  Changed: Fixed the computation of the finite element degree in
  VectorTools::project_boundary_values_curl_conforming_l2 to be based on the
  degree of the finite element selected for the application of the boundary
  values and not of the whole finite element.
  <br>
  (Winnifried Wollner, 2019/11/02)
 </li>

 <li>
  New: DoFTools::make_periodicity_constraints() can now be used to
  implement Bloch periodic conditions.
  <br>
  (Daniel Garcia-Sanchez, 2019/11/01)
 </li>

 <li>
  Changed: parallel::TriangulationBase::copy_triangulation does not
  copy the MPI communicator from the given source triangulation.
  <br>
  (Alexander Grayver, 2019/10/31)
 </li>

 <li>
  New: Member function hp::FECollection::max_degree() that returns the
  maximal polynomial degree of all finite elements in the collection.
  <br>
  (Marc Fehling, 2019/10/30)
 </li>

 <li>
  New: Added convenience functions CellAccessor::face_iterator_to_index()
  and TriaAccessor::child_iterator_to_index() for computing face and child numbers
  from iterators.
  <br>
  (David Wells, 2019/10/27)
 </li>

 <li>
  New: Added two new overloads of FEValues::reinit() that take face and subface
  iterators instead of face and subface numbers.
  <br>
  (David Wells, 2019/10/27)
 </li>

 <li>
  New: The Utilities::MPI::DuplicatedCommunicator helper class allows
  safe duplication and freeing of an MPI communicator.
  <br>
  (Timo Heister, 2019/10/22)
 </li>

 <li>
  Improved: VectorTools::interpolate() now works for a parallel::distributed::Triangulation
  <br>
  (Pasquale Claudio Africa, Doug Shi-Dong, 2019/10/24)
 </li>

 <li>
  New: new function Utilities::truncate_to_n_digits() that can be used to
  cutoff floating point numbers after a specified number of digits of
  accuracy (decimal places) in scientific floating point notation.
  <br>
  (Niklas Fehn, 2019/10/18)
 </li>

 <li>
  Refactor: improve implementation of Utilities::needed_digits() by using
  functionality from the standard library.
  <br>
  (Niklas Fehn, 2019/10/18)
 </li>

 <li>
  Changed: SolverControl returns number of iterations of ArpackSolver after solving
  for Eigenvalues. It also checks if the requested eigenvalues are supported by the
  solver, depending on symmetry.
  <br>
  (Roland Richter, 2019/10/17)
 </li>

 <li>
  New: Extend function GridGenerator::torus<3,3>() so that one can generate
  an open torus with an angle of 0 < phi <= 2*pi.
  <br>
  (Niklas Fehn, 2019/10/16)
 </li>

 <li>
  New: Provide cell iterators to interior and exterior cells of
  a macro face through the function MatrixFree::get_face_iterator().
  <br>
  (Martin Kronbichler, Peter Munch, Michal Wichrowski, 2019/10/16)
 </li>

 <li>
  New: The function DataOutInterface::write_vtu_with_pvtu_record() combines
  write_vtu() and write_pvtu_record() into a single function and generates
  the file extensions with counter, processor ID, and .vtu/.pvtu ending
  automatically.
  <br>
  (Niklas Fehn, 2019/10/14)
 </li>

 <li>
  New: The class VectorizedArray now provides STL-like member functions:
  VectorizedArray::size() returns the number of array elements and
  VectorizedArray:begin() as well as VectorizedArray::end() provide iterators
  enabling range-based iterations.
  <br>
  (Martin Kronbichler, Peter Munch, Daniel Arndt, 2019/10/12)
 </li>

 <li>
  Improved: Make the output of cells, (all non-internal) faces and (all non-internal)
  co-faces optional in GridOut::write_vtk.
  <br>
  (Peter Munch, Niklas Fehn, 2019/10/11)
 </li>

 <li>
  Fixed: Take in account p::d::Triangulation when using
  DoFHandler<dim, spacedim>::n_boundary_dofs() and its
  hp version. It now returns locally owned dofs instead
  of attempting a global count.
  <br>
  (Doug Shi-Dong, 2019/10/11)
 </li>

 <li>
  New: The function GridTools::compute_maximum_aspect_ratio() computes the
  maximum aspect ratio of a mesh with arbitrarily deformed elements.
  <br>
  (Niklas Fehn, 2019/10/08)
 </li>

 <li>
  Fixed: Make sure that merge_triangulations() does
  not forget about where the boundary is (the problem has
  actually been in create_triangulation(), which assigned
  numbers::internal_face_boundary_id to boundary faces)
  <br>
  (Sebastian Stark, 2019/10/08)
 </li>

 <li>
  Improved: CUDAWrappers::MatrixFree::cell_loop() can now detect whether the
  given vectors have the same parallel MPI partitioner as the matrix-free
  class, which circumvents vector copies and speeds up execution. This
  approach has also been applied to step-64.
  <br>
  (Martin Kronbichler, Peter Munch, 2019/09/27)
 </li>

 <li>
  Fixed: Manifold::get_new_point() previously ignored points weighted negatively.
  The original intent was to avoid to divide by a 0 weight, but the condition
  weight < tolerance resulted in ignoring negative weights.
  Additionally, the same function previously avoided get_intermediate_point()
  when the weight of of the first point was 0, but forgot to assign the new
  intermediate point to be the second point.
  The issues would mainly occur for derived Manifold classes that only implement
  the get_intermediate_point() function. Both issues are now fixed.
  <br>
  (Doug Shi-Dong, 2019/09/26)
 </li>

 <li>
  Improved: Evaluate/integrate the gradients in CUDAWrappers::EvaluatorTensorProduct
  with collocation method if both values and gradients are requested.
  <br>
  (Peter Munch, 2019/09/24)
 </li>

 <li>
  New: Provide new methods in CUDAWrappers::FEEvaluation that do neither take the
  local dof index nor the quadrature point but recompute them based on the
  thread id.
  <br>
  (Peter Munch, 2019/09/24)
 </li>

 <li>
  Fixed: When called with update_volume_elements, MappingFEField had an
  assertion that accidentally checked the wrong field's size. This is
  now fixed.
  <br>
  (Doug Shi-Dong, 2019/09/23)
 </li>

 <li>
  Improved: Unify compute_vertices_with_ghost_neighbors in
  parallel::TriangulationBase, parallel::distributed::Triangulation,
  and parallel::fullydistributed::Triangulation..
  <br>
  (Peter Munch, 2019/09/21
 </li>

 <li>
  New: The class DiscreteTime is created to keep track of time
  during a time-dependent simulation.
  <br>
  (Reza Rastak, 2019/09/20)
 </li>

 <li>
  Changed: The storage type of MGLevelGlobalTransfer::copy_indices,
  MGLevelGlobalTransfer::copy_indices_global_mine,
  MGLevelGlobalTransfer::copy_indices_level_mine has been changed to a
  std::vector of Table<2,unsigned int>, in order to simplify some copy
  operations.
  <br>
  (Martin Kronbichler, 2019/09/20)
 </li>

 <li>
  Improved: MGTools::make_boundary_list() has been made faster.
  <br>
  (Martin Kronbichler, 2019/09/19)
 </li>

 <li>
  Improved: IndexSet::add_indices() could previously have quadratic complexity
  in the number of ranges if many ranges were added. This is now fixed. Note,
  however, that calling IndexSet::add_indices() many times still is of quadratic
  complexity.
  <br>
  (Martin Kronbichler, 2019/09/19)
 </li>

 <li>
  New: MGTransferMatrixFree::build() now takes an optional second argument to
  specify partitioner objects for LinearAlgebra::distributed::Vector types in
  use elsewhere. If the ghost indices required by the multigrid algorithm are
  contained in the ghosts of the given partitioners, they can be used, saving
  some vector copy operations during MGTransferMatrixFree::restrict_and_add()
  and MGTransferMatrixFree::prolongate().
  <br>
  (Martin Kronbichler, 2019/09/18)
 </li>

 <li>
  Improved: The ghost setup of Utilities::MPI::Partitioner has been rewritten
  using Utilities::MPI::ConsensusAlgorithm to avoid global all-to-all
  communication, speeding up computations on many MPI ranks.
  <br>
  (Peter Munch, Martin Kronbichler, 2019/09/17)
 </li>

 <li>
  Improved: The gzip compression previously applied by
  GridTools::exchange_cell_data_to_ghosts() has been disabled because network
  throughput is typically higher than the rate by which gzip is able to
  compress. This speeds up DoFHandler::distribute_dofs(), among others.
  <br>
  (Martin Kronbichler, 2019/09/16)
 </li>

 <li>
  Added: Add a new particle generation function that can generate particles at
  random positions controlled by a probability density function.
  <br>
  (Rene Gassmoeller, 2019/09/13)
 </li>

 <li>
  New: Allow to create TriangulationDescription::Description
  by groups, in which only one process actually creates
  TriangulationDescription::Descriptions and distributes them within the group.
  <br>
  (Peter Munch, 2019/09/13)
 </li>

 <li>
  New: Add Boost serialization function to dealii::TriangulationDescription::CellData and
  dealii::TriangulationDescription::Description.
  <br>
  (Peter Munch, 2019/09/13)
 </li>

 <li>
  New: Add Boost serialization function to CellData.
  <br>
  (Peter Munch, 2019/09/13)
 </li>

 <li>
  Fixed: DoFHandler::distribute_mg_dofs() for distributed triangulations has
  been made up to ten times faster by replacing an inefficient function for the
  exchange of ghosted level degrees of freedom.
  <br>
  (Martin Kronbichler, 2019/09/13)
 </li>

 <li>
  New: MGConstrainedDoFs now allows the user to define constraints
  on each multigrid level to be applied during prolongation in a
  matrix-free transfer.
  <br>
  (Conrad Clevenger, 2019/09/13)
 </li>

 <li>
  New: Add method is_parent_of to CellId.
  <br>
  (Peter Munch, 2019/09/11)
 </li>

 <li>
  New: There is a new variant of MatrixFree::cell_loop() which takes two
  std::function objects with ranges on the locally owned degrees of freedom, one
  with work to be scheduled before the cell operation first touches some
  unknowns, and one after the cell operation last touches them. The goal of
  these functions is to bring vector operations close to the time when the
  vector entries are accessed by the cell operation which increases the cache
  hit rate of modern processors (temporal locality).
  <br>
  (Martin Kronbichler, 2019/09/10)
 </li>

 <li>
  New: compare_and_apply_mask() allows to use a ternary operator idiom with
  VectorizedArray that is based on generating a mask via component-wise
  comparison. Such a computational idiom is useful as an alternative to
  branching whenever the control flow itself would depend on (computed) data
  (which is not possible for vectorized data).
  <br>
  (Matthias Maier, 2019/09/05)
 </li>

 <li>
  Fixed: Implementing LinearAlgebra::Vector::has_ghost_elements() improves compatibility
  between vector types.
  <br>
  (Daniel Arndt, 2019/09/05)
 </li>

 <li>
  Added: GeometryInfo<dim>::unit_normal_vector provides outward-facing unit
  vectors for each face of the reference cell in the form of Tensor<1, dim>.
  <br>
  (Reza Rastak, 2019/09/04)
 </li>

 <li>
  Fixed: PreconditionChebyshev can be used with LinearAlgebra::distributed::Vector
  with MemorySpace::CUDA. In particular, PreconditionChebyshev can set constrained
  indices to zero in the initial guess.
  <br>
  (Daniel Arndt, 2019/09/03)
 </li>

 <li>
  New: Add a function, which creates a TriangulationDescription::Description from a
  serial and a distributed triangulation. This data structure can be used to initialize
  parallel::fullydistributed::Triangulation.
  <br>
  (Peter Munch, 2019/09/02)
 </li>

 <li>
  Bugfix: In PolynomialsBDM, the degree() function is corrected to return a
  polynomial degree of k+1, instead of its previous value of k.
  <br>
  (Graham Harper, 2019/09/02)
 </li>

 <li>
  Fixed: Reading the variable <code>debug</code> using ParameterHandler for the class
  Algorithms::Newton<VectorType> is fixed within the method
  Algorithms::Newton<VectorType>::parse_parameters().
  <br>
  (Reza Rastak, 2019/08/31)
 </li>

 <li>
  Fixed: CUDAWrappers::MatrixFree::cell_loop() now works on adapted mesh even in
  two dimensions with finite elements of even degree.
  <br>
  (Bruno Turcksin, 2019/08/31)
 </li>

 <li>
  New: Added a GridGenerator::generate_from_name_and_arguments() that leverages the MutableBind class
  to parse two strings, and generate a grid from them interpreting one as a GridGenerator function
  name to call, and the other as a tuple of arguments to pass to the function.
  <br>
  (Luca Heltai, 2019/08/27)
 </li>

 <li>
  Changed: No of lines reported in step-3 documentation results section is corrected to two
  <br>
  (Krishnakumar Gopalakrishnan, 2019/08/27)
 </li>

 <li>
  Changed: In Step-3 docs, clarify wording on which objects belong to SparseMatrix & Vector classes.
  <br>
  (Krishnakumar Gopalakrishnan, 2019/08/27)
 </li>

 <li>
  Bugfix: For parallel::distributed::Triangulation objects, p4est determines
  which cells will be refined and coarsened. Thus, refinement flags are not
  a good measure to predict refinement behavior for transferring active finite
  element indices. We now provide a matching coarsening strategy to the
  parallel::distributed::CellDataTransfer object responsible for their transfer.
  This fixes assertions being triggered in
  parallel::distributed::SolutionTransfer::interpolate().
  <br>
  (Marc Fehling, 2019/08/27)
 </li>

 <li>
  New: There is a new constructor for DiagonalMatrix that immediately initializes
  the underlying vector.
  <br>
  (Daniel Arndt, 2019/08/23)
 </li>

 <li>
  Changed: step-3 now uses VTK as the output file format.
  <br>
  (Wolfgang Bangerth, 2019/08/23)
 </li>

 <li>
  New: method try_get_data() added to CellDataStorage in both const and
  non-const forms. It returns a std_cxx17::optional which determines whether a
  cell has a data associated with it.
  <br>
  (Reza Rastak, 2019/08/21)
 </li>

 <li>
  Fixed: OpenCASCADE::interpolation_curve() returned a TopoDS_Shape that was not properly closed. This is
  now fixed.
  <br>
  (Luca Heltai, 2019/08/21)
 </li>

 <li>
  Fixed: Make sure that begin() and end() of PETScWrappers::MatrixBase can only be called on a processor owning all rows of the matrix. Also fix a bug which caused begin() to fail when the first row(s) of the matrix is(are) empty.
  <br>
  (Sebastian Stark, 2019/08/19)
 </li>

 <li>
  New: Implement DoFAccessor::set_mg_dof_indices and get_mg_dof_indices for 1D.
  <br>
  (Peter Munch, 2019/08/18)
 </li>

 <li>
  New: Introduce new class parallel::DistributedTriangulationBase
  between parallel::TriangulationBase and parallel::distributed::Triangulation.
  <br>
  (Peter Munch, 2019/08/17)
 </li>

 <li>
  Fixed: Forward declarations in header files are prevented from confusing
  doxygen.
  <br>
  (Reza Rastak, 2019/08/16)
 </li>

 <li>
  New: Enable internal::DoFHandlerImplementation::Policy::ParallelDistributed for 1D.
  <br>
  (Peter Munch, 2019/08/16)
 </li>

 <li>
  New: Add method to get the coarse-grid cell from CellID.
  <br>
  (Peter Munch, 2019/08/15)
 </li>

 <li>
  Added: Add a Particles::Generators namespace and a first implementation that
  generates particles at regular position in the reference domain.
  <br>
  (Rene Gassmoeller, 2019/08/14)
 </li>

 <li>
  New: Add CellAccessor::face_iterators() and DoFCellAccessor::face_iterators()
  which return an iterator range over the faces of a cell.
  <br>
  (Peter Munch, Wolfgang Bangerth, Daniel Arndt, 2019/08/14)
 </li>

 <li>
  New: Generalize internal::DoFHandlerImplementation::Policy::ParallelDistributed such that it uses
  the new definition of CellId.
  <br>
  (Peter Munch, 2019/08/14)
 </li>

 <li>
  New: The definition of CellID has been modified so that it depends on the
  unique coarse-cell id and no longer on the coarse-cell index.
  <br>
  (Peter Munch, 2019/08/13)
 </li>

 <li>
  Fixed: MatrixFree::loop() and MatrixFree::cell_loop() forgot to exchange
  ghosted vector entries located on the first active cell of an MPI rank for
  continuous elements when face integrals are activated. This is now fixed.
  <br>
  (Martin Kronbichler, Fabian Castelli, 2019/08/11)
 </li>

 <li>
  Improved: Make a large number of prose, grammar, and various other improvements
  to the library.
  <br>
  (Ashna Aggarwal, Mario Zepeda Aguilar, Manaswinee Bezbaruah, Bruno Blais, Kirana Bergstrom, Katherine Cosburn, Luel Emishaw, Rebecca Fildes, Andres Galindo, Brandon Gleeson, Bang He, Nicole Hayes, Jiuhua Hu, Lise-Marie Imbert-Gerard, Manu Jayadharan, Marie Kajan, Charu Lata, Adriana Morales Miranda, Emily Novak, Rebecca Pereira, Geneva Porter, Irabiel Romero, Tonatiuh Sanchez-Vizuet, Sara Tro, 2019/08/08)
 </li>

 <li>
  New: Added the function GridGenerator::subdivided_hyper_L().
  <br>
  (Mae Markowski, 2019/08/10)
 </li>

 <li>
  Improved: parallel::distributed::ContinuousQuadratureDataTransfer allows the
  number of data values per quadrature points vary across the mesh. It also allows
  a cell to not have any associated quadrature data during refinement.
  <br>
  (Reza Rastak, Marc Fehling, Peter Munch, 2019/08/09)
 </li>

 <li>
  New: Add a possible extension to step-2.
  <br>
  (Wenjuan Zhang, 2019/08/08)
 </li>

 <li>
  Added: GridGenerator::eccentric_hyper_shell() a volume contained between two circles or spheres centered around two different points.
  <br>
  (Melanie Gerault, 2019/08/08)
 </li>

 <li>
  Fixed: Fixed two functions in GridTools::Cache class, namely
  GridTools::Cache::get_used_vertices_rtree() and
  GridTools::Cache::get_cell_bounding_boxes_rtree() where the update flags were
  not cleared. Fixed that by adding a line to both functions at the end of the if
  clause checking the flag.
  <br>
  (Manu Jayadharan, 2019/08/08)
 </li>

 <li>
  New: Added Utilities::MutableBind class and Utilities::mutable_bind() method.
  <br>
  (Luca Heltai, Matthias Maier, 2019/08/08)
 </li>

 <li>
  New: Added std_cxx17::apply() method.
  <br>
  (Luca Heltai, 2019/08/08)
 </li>

 <li>
  New: Add a new parameter to CUDAWrappers::MatrixFree::AdditionalData to choose
  between graph coloring and atomics in CUDAWrappers::MatrixFree::cell_loop().
  The atomic implementation is faster on Pascal and newer architectures.
  <br>
  (Bruno Turcksin, 2019/08/08)
 </li>

 <li>
  Fixed: Bug in CellDataStorage that produced error when applying adaptive
  refinement in multi-material grids.
  <br>
  (Reza Rastak, Peter Munch, 2019/08/07)
 </li>

 <li>
  Changed: The saturation equation in step-21 is explained further for
  consistency with literatures in porous media transport.
  <br>
  (Omotayo Omosebi, 2019/08/07)
 </li>

 <li>
  Improved: CellDataStorage::get_data() now checks whether we are retrieving
  using a data type that matches the type we used to initialize the cell data.
  <br>
  (Reza Rastak, 2019/08/06)
 </li>

 <li>
  New: Added BlockMatrixBase::frobenius_norm().
  (Jonathan Robey, 2019/08/06)
 </li>

 <li>
  Improved: Unified in-line journal citations in steps 1, 6, and 18; Bibtex
  style formatting in a central bibliography.
  <br>
  (Adam Lee, 2019/08/06)
 </li>

 <li>
  Changed: The step-5 program no longer uses EPS to output the solutions,
  but the far more modern VTK/VTU format.
  <br>
  (Wolfgang Bangerth, 2019/08/03)
 </li>

 <li>
  Changed: The step-1 program no longer uses EPS to output the meshes,
  but the far more modern SVG format.
  <br>
  (Wolfgang Bangerth, 2019/08/03)
 </li>

 <li>
  Improved: DoFTools::make_periodicity_constraints() can now be
  used with std::complex<double> numbers.
  <br>
  (Daniel Garcia-Sanchez, 2019/08/02)
 </li>

 <li>
  New: Amend VectorizedArray with a default constructor and
  a constructor taking a scalar value.
  <br>
  (Peter Munch, 2019/07/29)
 </li>

 <li>
  New: The content of a VectorizedArray object can be printed
  using the output operator.
  <br>
  (Daniel Arndt, 2019/07/29)
 </li>

 <li>
  Fixed: AffineConstraints::set_zero() can be used with
  LinearAlgebra::distributed::Vector with MemorySpace::CUDA.
  <br>
  (Daniel Arndt, 2019/07/28)
 </li>

 <li>
  New: Enable periodicity for 1D triangulations (only for dealii::DoFHandler).
  <br>
  (Peter Munch, 2019/07/27)
 </li>

 <li>
  New: The function LinearAlgebra::distributed::Vector::import() allows to
  copy the data from a LinearAlgebra::distributed::Vector of a given MemorySpace
  to another LinearAlgebra::distributed::Vector of a different MemorySpace.
  <br>
  (Bruno Turcksin, 2019/07/24)
 </li>

 <li>
  Improved: Add assignment operator for two Point objects with different
  underlying scalar types.
  <br>
  (Daniel Garcia-Sanchez, 2019/07/23)
 </li>

 <li>
  Fixed: When compiling with CUDA support, the call to
  LinearAlgebra::distributed::Vector::compress() would failed if the vector's
  memory space was MemorySpace::Host, the header
  deal.II/lac/la_parallel_vector.templates.h was included, and the call was done
  in a .cu file.
  <br>
  (Bruno Turcksin, 2019/07/23)
 </li>

 <li>
  New: MappingFEField can now also be initialized on multigrid levels by handing
  in an appropriately sized vector of position vectors for the levels.
  <br>
  (Martin Kronbichler, Johannes Heinz, 2019/07/22)
 </li>

 <li>
  New: CUDAWrappers::MatrixFree::initialize_dof_vector() can be used to
  initialize vectors suitable for the parallel partitioning of the degrees of
  freedoms similar to MatrixFree::initialize_dof_vector().
  <br>
  (Daniel Arndt, 2019/07/19)
 </li>

 <li>
  Improved: constexpr evaluation is enabled in SymmetricTensor and TableIndices.
  <br>
  (Reza Rastak, 2019/07/13)
 </li>

 <li>
  New: Provide a list of all possible VectorizedArray<Number, width> instances for a given
  hardware and optimization level. The list can be used during template instantiation.
  <br>
  (Peter Munch, 2019/07/13)
 </li>

 <li>
  Improved: The additional roots of the HermiteLikeInterpolation with degree $p$
  greater than four have been switched to the roots of the Jacobi polynomial
  $P^{(4,4)}_{p-3}$, making the interior bubble functions $L_2$ orthogonal and
  improving the conditioning of interpolation slightly.
  <br>
  (Martin Kronbichler, 2019/07/12)
 </li>

 <li>
  Improved: FE_PolyTensor now supports finite elements with different mapping
  types for each shape function. The member variable mapping_type is now of
  type std::vector<MappingType>. If the vector contains only one MappingType,
  then the same mapping is applied to every shape function. If a list of
  MappingTypes is provided, mapping_type[i] will be applied to shape function i.
  <br>
  (Graham Harper, 2019/07/12)
 </li>

 <li>
  New: It is possible to access the value type Number of VectorizedArray<Number,width>
  directly via VectorizedArray<Number, width>::value_type.
  <br>
  (Peter Munch, 2019/07/10)
 </li>

 <li>
  New: CUDAWrappers::FEEvaluation::get_dof_values() and
  CUDAWrappers::FEEvaluation::submit_dof_values() provide access to the values
  stored for the degrees of freedom.
  <br>
  (Daniel Arndt, 2019/07/09)
 </li>

 <li>
  New: The library contains a new manual page describing how to package deal.II
  for distribution in a Unix-like distribution.
  <br>
  (David Wells, 2019/07/07)
 </li>

 <li>
  New: An additional optional template parameter VectorizedArrayType has been
  added to the MatrixFree framework. It lets the user switch between no
  vectorization (only auto vectorization) and vectorization over elements.
  <br>
  (Peter Munch, 2019/07/06)
 </li>

 <li>
  Improved: All member functions and free functions related to Point
  that can be used in CUDA code so far have been annotated accordingly.
  <br>
  (Daniel Arndt, 2019/07/05)
 </li>

 <li>
  New: A second template parameter has been added to VectorizedArray. This template
  parameter optionally controls the number of vectorization lanes to be used, i.e.,
  which instruction set architecture extension to be used, at the time the application
  is built.
  <br>
  (Peter Munch, 2019/07/04)
 </li>

 <li>
  Improved: All member functions and free functions related to Tensor
  that can be used in CUDA code so far have been annotated accordingly.
  <br>
  (Daniel Arndt, 2019/07/03)
 </li>

 <li>
  Fixed: CUDAWrappers::MatrixFree::cell_loop() would deadlock on adaptively
  refined meshes when using Volta GPU.
  <br>
  (Bruno Turcksin, 2019/07/02)
 </li>

 <li>
  New: Class CellDataTransfer has been introduced to transfer cell by cell
  data across meshes in case of refinement/serialization.
  <br>
  (Marc Fehling, 2019/06/26)
 </li>

 <li>
  Improved: The functions vectorized_load_and_transpose() and
  vectorized_transpose_and_store() for AVX-512 were rewritten to better utilize
  execution units with two load units and one shuffle/swizzle unit as common on
  recent Intel CPUs.
  <br>
  (Martin Kronbichler, 2019/06/24)
 </li>

 <li>
  Improved: The function project_boundary_values_curl_conforming_l2() can now be
  used with std::complex<double> numbers.
  <br>
  (Daniel Garcia-Sanchez, 2019/06/22)
 </li>

 <li>
  New: The function schur_product was added for the Tensor class,
  allowing the entrywise multiplication of tensor objects of general rank.
  <br>
  (Roland Richter, 2019/06/10)
 </li>

 <li>
  New: Namespace hp::Refinement offering decision tools for p adaptivity.
  <br>
  (Marc Fehling, 2019/06/07)
 </li>

 <li>
  Unified DoFTools::extract_dofs() for DoFHandler and hp::DoFHandler.
  <br>
  (Daniel Arndt, Mathias Anselmann, 2019/06/07)
 </li>

 <li>
  New: The class MappingQCache implements a cache for the points generated by
  MappingQGeneric and derived classes to speed up consumers of the mapping for
  expensive manifolds.
  <br>
  (Martin Kronbichler, 2019/06/06)
 </li>

 <li>
  New: The function Utilities::MPI::compute_index_owner() provides the owner of
  the ghost degrees of freedom relative to a set of locally owned ones.
  <br>
  (Peter Munch, 2019/06/03)
 </li>

 <li>
  Improved: A number of places in deal.II previously called `MPI_Allgather`
  followed by some local sums to compute prefix sums. Now the more appropriate
  `MPI_Exscan` is used instead, which improves scalability of setup routines on
  more than 10k processors.
  <br>
  (Martin Kronbichler, Peter Munch, 2019/06/03)
 </li>

 <li>
  New: Add functions SparseMatrix::print_as_numpy_arrays() and
  LinearAlgebra::Vector::print_as_numpy_array() to output the
  object in a format readable with numpy.loadtxt.
  <br>
  (Bruno Turcksin, 2019/06/02)
 </li>

 <li>
  New: Tensor<rank, dim, Number> can be constructed and manipulated in constexpr
  setting (if supported by the compiler).
  <br>
  (Reza Rastak, 2019/05/26)
 </li>

 <li>
  Fixed: FullMatrix can now be packed into std::vector and unpacked from it using
  std::copy and the begin() and end() iterators.
  <br>
  (Reza Rastak, 2019/05/24)
 </li>

 <li>
  Fixed: TrilinosWrappers::SparsityPattern can deal with empty column maps.
  <br>
  (Daniel Arndt, Mathias Anselmann, 2019/05/22)
 </li>

 <li>
  Improved: The function GridGenerator::torus() run in 3D (volume mesh) would
  previously only use a single cell to represent the poloidal shape, which leads
  to singular mappings similar to how a circle degenerates when meshed with a
  single cell. The poloidal shape is now represented with 5 cells just like the
  circle. Furthermore, the function GridGenerator::torus() has gained an
  optional argument to control the number of cell layers in the toroidal
  direction. The default manifold set to the torus has also been improved: Now,
  the TorusManifold is applied on the surface, a CylindricalManifold to the
  middle cells in toroidal coordinates, and a TransfiniteInterpolationManifold
  on the cells between the surface and the inner ring. This leads to an
  excellent mesh quality for all supported scenarios.
  <br>
  (Niklas Fehn, Martin Kronbichler, 2019/05/16)
 </li>

 <li>
  Fixed: Fix a bug in tutorial step-46, which caused the y-displacement (z-displacement in 3d) to be slightly unsymmetric and rectify the interface traction term.
  <br>
  (Sebastian Stark, 2019/05/15)
 </li>

 <li>
  New: Two new methods NonMatching::compute_coupling_sparsity_pattern() and
  NonMatching::compute_coupling_mass_matrix() allow to construct the coupling between arbitrary grids,
  using a convolution kernel.
  <br>
  (Luca Heltai, Wenyu Lei, 2019/04/19)
 </li>

 <li>
  New: Add FESeries::Legendre::get_n_coefficients_per_direction() and
  FESeries::Fourier::get_n_coefficients_per_direction() to retrieve
  the number of coefficients in each direction. Also add an Assert in
  FESeries::Legendre::calculate() and FESeries::Fourier::calculate() to
  check the dimension of the table to store coefficients.
  <br>
  (Denis Davydov, 2018/12/27)
 </li>

 <li>
  New: The InterpolatedUniformGridData::gradient() function was not
  previously implemented. It is now.
  <br>
  (Bob Myhill, Anne Glerum, Wolfgang Bangerth, 2018/09/27)
 </li>

</ol>

*/
