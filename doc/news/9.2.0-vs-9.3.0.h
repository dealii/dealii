// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
@page changes_between_9_2_0_and_9_3_0 Changes between Version 9.2.0 and 9.3.0

<p>
This is the list of changes made between the release of deal.II version
9.2.0 and that of 9.3.0. All entries are signed with the names of the
author.
</p>
<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="920-930-incompatible"></a>
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
  Removed: An additional overload of GridTools::cell_measure(), which only
  existed to aid generic programming and only threw an exception, has been
  removed. If you use this function then the best fix is to ensure that
  you are not explicitly providing a template parameter to the function itself.
  <br>
  (David Wells, 2021/05/16)
 </li>

 <li>
  Renamed: The function ParticleHandler::locally_relevant_ids() has been deprecated.
  Please use the new function ParticleHandler::locally_owned_particle_ids() instead.
  <br>
  (Peter Munch, 2021/05/01)
 </li>

 <li>
  Deprecated: The GridReordering class as well as
  Triangulation::create_triangulation_compatibility have been deprecated.
  These functions use the old-style (before 5.2) numbering and have been
  unofficially deprecated since 2005.
  <br>
  (David Wells, 2021/04/22)
 </li>

 <li>
  Removed: The class parallel::distributed::ErrorPredictor has been
  removed. Use the function hp::Refinement::predict_error() in combination
  with parallel::distributed::CellDataTransfer instead. Please consult the
  documentation of hp::Refinement::predict_error() for instructions.
  <br>
  (Marc Fehling, 2021/04/21)
 </li>

 <li>
  Changed: The order in which vertices are read from GridIn::read_unv() has
  changed. As a result the vertex numbers on each cell are now slightly different.
  <br>
  (David Wells, 2021/04/13)
 </li>

 <li>
  Changed: The order of parameters has been switched in MappingQCache::initialize(),
  with the Mapping now being the first argument.
  <br>
  (Peter Munch, 2020/05/23)
 </li>

 <li>
  Changed: The various parse_arguments(), to_string(), and to_value() functions
  related to the Pattern namespace take their last argument by `const&` instead of
  `const std::unique_ptr<>&`.
  <br>
  (Daniel Arndt, 2021/02/22)
 </li>

 <li>
  Fixed: GridOut::write_vtk output the vertices of hexahedra in the wrong order,
  which resulted in visualization programs computing all cell volumes as negative.
  Fixing this required that we significantly update the way VTK input is read in
  to deal.II as well (so that GridOut composed with GridIn remains idempotent). As
  a result, the order of cells, faces, and edges is now different for meshes read
  from VTK files.
  <br>
  (David Wells, 2021/02/20)
 </li>

 <li>
  Changed: The FE_Q_Base class took a description of the polynomial
  space as a template argument. But that is not necessary: It is
  entirely sufficient to pass this information to the constructor in the
  form of a regular constructor argument. This has been changed, and the
  class has therefore lost its first template argument.
  <br>
  (Wolfgang Bangerth, 2021/02/03)
 </li>

 <li>
  Deprecated: ARKode does no longer need a `reinit_vector` function.
  <br>
  (Sebastian Proell, 2021/02/01)
 </li>

 <li>
  Added a simple mesh for debugging purposes that consists
  of two cubes where one of them can be chosen to have one face
  either flipped or rotated.
  <br>
  (Konrad Simon, 2021/01/28)
 </li>

 <li>
  Updated: deal.II dropped support for Sundials older than version 3.0.0.
  <br>
  (Matthias Maier, Luca Heltai, 2021/01/25)
 </li>

 <li>
  Updated: deal.II now requires CMake version 3.1.0 or newer to configure
  <br>
  (Matthias Maier, 2021/01/25)
 </li>

 <li>
  Deprecated: All `LinearAlgebra::distributed::Vector::zero_out_ghosts()` has been deprecated.
  Use `LinearAlgebra::distributed::Vector::zero_out_ghost_values()` instead.
  <br>
  (Peter Munch, 2021/01/04)
 </li>

 <li>
  Patch to fix issue#7970 for 'RaviartThomas<3>(degree)' elements on meshes
  that contain cells with faces that are flipped and/or not in standard
  orientation. The patch provides a possible way to fix other FE classes inheriting
  from 'FE_PolyTensor' but does not implement it yet. It does not change the interface
  to FE classes.
  <br>
  (Konrad Simon, 2020/12/22)
 </li>

 <li>
  Changed: The interface of FE::fill_fe_face_values() now accepts instead
  of a Quadrature instance a hp::QCollection instance, enabling the evaluation of
  shape functions for different quadrature formulas on different faces for
  FiniteElement classes with supports for this feature (e.g. FE_Q, FE_SimplexP).
  <br>
  (Peter Munch, 2020/12/12)
 </li>

 <li>
  Deprecated: The function CellId::to_cell() has been replaced by
  Triangulation::create_cell_iterator().
  <br>
  (Marc Fehling, 2020/12/12)
 </li>

 <li>
  Changed: GridTools::find_active_cell_around_point() no longer throws an exception, but returns an invalid iterator.
  User codes should now check that the returned cell is valid instead of relying on exceptions, instead of trying
  to catch exceptions that this function may have thrown.
  <br>
  (Luca Heltai, 2020/11/30)
 </li>

 <li>
  Removed: The file `affine_constraints.h` had a class called
  `IsBlockMatrix` that, however, also allowed to check whether its
  template argument is in fact a block sparsity pattern. The class was
  likely never intended to be in the public interface, and has been
  moved to internal::AffineConstraints::IsBlockMatrix and now also
  internal::AffineConstraints::IsBlockSparsityPattern.
  <br>
  (Wolfgang Bangerth, 2020/11/29)
 </li>

 <li>
  Deprecated: The DoFHandlerType template argument for the DataOutStack
  class has been deprecated. Use DataOutStack<dim, spacedim> instead.
  <br>
  (Marc Fehling, 2020/11/20)
 </li>

 <li>
  Deprecation announcement: The template arguments of the following
  classes will change in a future release:
  <ul>
    <li>SolutionTransfer
    <li>parallel::distributed::SolutionTransfer
    <li>Functions::FEFieldFunction
    <li>DataOut
    <li>DataOut_DoFData
    <li>DataOutFaces
    <li>DataOutRotation
  </ul>
  If for some reason, you need a code that is compatible with deal.II
  9.3 and the subsequent release, a Legacy namespace has been introduced
  with aliases of these classes with their original interface. You can
  make the following substitutions to your code for each of the affected
  classes:
  <ul>
    <li>X &rarr; Legacy::X
  </ul>
  To perform this substitution automatically, you may use a *search and
  replace* script like the following made for the *bash* shell:
  @code{.sh}
  classes=(SolutionTransfer parallel::distributed::SolutionTransfer Functions::FEFieldFunction DataOut DataOut_DoFData DataOutFaces DataOutRotation)
  for c in \${classes[@]}; do
    find /path/to/your/code -type f -exec sed -i -E "/(\w\${c}|\${c}[^<]|Particles::\${c}|distributed::\${c}|^\s*(\/\/|\*))/! s/\${c}/Legacy::\${c}/g" {} \;
  done
  @endcode
  (Marc Fehling, 2020/11/20)
 </li>

 <li>
  Deprecated: The DoFHandlerType template argument for the functions
  DataPostprocessorInputs::CommonInputs::set_cell() and
  DataPostprocessorInputs::CommonInputs::get_cell() has been deprecated.
  These functions now only work with the basic `DoFHandler` class. As a
  consequence, the get_cell() function requires an additional dim
  template. For example, write: cell=common_inputs.template get_cell<dim>()
  <br>
  (Marc Fehling, 2020/11/15)
 </li>

 <li>
  Deprecated: The DoFHandlerType template argument for the MappingFEField
  class has been deprecated. Use MappingFEField<dim, spacedim, VectorType>
  instead.
  <br>
  (Marc Fehling, 2020/11/15)
 </li>

 <li>
  Deprecated: The variant of the function MatrixFree::get_dof_handler
  expecting a DoFHandlerType template has been deprecated. Use the
  template-less variant returning a DoFHandler instead.
  <br>
  (Marc Fehling, 2020/11/13)
 </li>

 <li>
  Deprecated: All `MatrixFree::reinit()` functions without `Mapping` as argument
  have been deprecated. Use the functions that explicit provide the Mapping instead.
  <br>
  (Peter Munch, 2020/11/11)
 </li>

 <li>
  Deprecated: The operator FiniteElement::operator[] has been deprecated.
  Use DoFHandler::get_fe() with a specified index instead of code like
  dof_handler->get_fe()[index].
  <br>
  (Marc Fehling, 2020/11/11)
 </li>

 <li>
  Deprecated: The initialization interface for the DoFHandler class has
  changed. DoFHandler::initialize() and DoFHandler::set_fe() have been
  deprecated. Instead, use the constructor DoFHandler::DoFHandler(tria) or
  DoFHandler::DoFHandler() followed by DoFHandler::reinit(tria) to attach
  a Triangulation tria, and DoFHandler::distribute_dofs() to enumerate
  degrees of freedom.
  <br>
  (Marc Fehling, 2020/10/30)
 </li>

 <li>
  Changed: MatrixFree::loop() only accepts LinearAlgebra::distributed::Vector arguments
  that have been initialized via MatrixFree::initialize_dof_vector()
  and are as a consequence globally compatible with the
  Utilities::MPI::Partitioner within internal::MatrixFreeFunctions::DoFInfo.
  <br>
  (Peter Munch, Martin Kronbichler, 2020/10/21)
 </li>

 <li>
  Changed: The return type of the method Mapping::get_vertices() has been changed
  from std::array to boost::container::small_vector.
  <br>
  (Peter Munch, 2020/10/16)
 </li>

 <li>
  Deprecated: The QTrapez class, poorly named because the proper English
  term is "trapezoidal quadrature rule", has been renamed to QTrapezoid,
  and the class with the old name has been deprecated.
  <br>
  (Wolfgang Bangerth, 2020/09/28)
 </li>

 <li>
  Replaced: Python wrapper for 'merge_triangulations' with a more generic equivalent.
  <br>
  (Alexander Grayver, 2020/09/01)
 </li>

 <li>
  Removed: CUDA 9, 10.0, and 10.1 are not supported anymore.
  <br>
  (Bruno Turcksin, 2020/08/05)
 </li>

 <li>
  Changed: The template arguments of the classes DoFAccessor and DoFCellAccessor have changed.
  The template argument DoFHandlerType has been replaced by dimension and
  space dimension.
  <br>
  (Peter Munch, 2020/06/24)
 </li>

 <li>
  Deprecated: The functions MatrixFree::reinit(), which take
  a vector of hp::DoFHandlers, have been deprecated. Users are asked
  to provide vectors of DoFhandlers, which may contain hp::DoFHandlers. This is
  possible now since hp::DoFHandler is deriving from DoFHandler.
  <br>
  (Peter Munch, 2020/06/03)
 </li>

 <li>
  Removed: The deprecated class MGCoarseGridLACIteration has been removed.
  <br>
  (Daniel Arndt, 2020/06/12)
 </li>

 <li>
  Removed: The header file `deal.II/grid/tria_object.h` has been
  removed. It was only used for internal purposes.
  <br>
  (Wolfgang Bangerth, 2020/06/05)
 </li>

 <li>
  Changed: The binary representation of the Triangulation and DoFHandler classes
  created by the function `save()` with release 9.2 cannot be read anymore with
  `load()` due to major internal changes
  of these classes. This change also affects, i.a., the functions
  `GridIn::read_vtu()` and `GridOut::write_vtu()`.
  <br>
  (Peter Munch, 2020/06/03)
 </li>

 <li>
  Removed: The deprecated bindings to the legacy NETCDF C++ library have been
  removed.
  <br>
  (David Wells, 2020/05/27)
 </li>

 <li>
  Removed: The deprecated bindings to nanoflann have been removed.
  <br>
  (David Wells, 2020/05/27)
 </li>

 <li>
  Changed: The ThreadLocalStorage class has been reimplemented with C++14 STL
  primitives and does not depend on the TBB library any more. With that the
  obscure ThreadLocalStorage::get_implementation() function that exposed the
  underlying TBB container has been removed.
  <br>
  (Matthias Maier, 2020/05/23)
 </li>

 <li>
  Removed: The Threads::Task class had an `operator==` that allowed
  comparing objects of this type for equality. This operator has been
  removed. If you want to store objects of this kind in a collection
  that requires this kind of operator (say, `std::set`), then you
  probably can't do so any more in a reasonable way. However, this is
  exactly what the Threads::TaskGroup class is there for.
  <br>
  (Wolfgang Bangerth, 2020/05/26)
 </li>

 <li>
  Deprecated: The class hp::DoFHandler has been deprecated, since the DoFHandler
  has been extended with its functionalities.
  <br>
  (Peter Munch, 2020/05/23)
 </li>

 <li>
  Removed: The following preprocessor definitions have been removed from
  config.h.in: DEAL_II_NOEXCEPT, DEAL_II_USE_MT_POSIX,
  DEAL_II_USE_MT_POSIX_NO_BARRIERS
  <br>
  (Matthias Maier, 2020/05/23)
 </li>

 <li>
  Removed: The deprecated classes Threads::Mutex::ScopedLock,
  Threads::ConditionVariable, and deprecated functions
  Threads::Mutex::acquire(), Threads::Mutex::release(),
  Threads::n_existing_threads(), Threads::this_thread_id() have been removed.
  <br>
  (Matthias Maier, 2020/05/22)
 </li>

 <li>
  Removed: The <code>DEAL_II_WITH_CXX14</code> and
  <code>DEAL_II_WITH_CXX17</code> configuration options have been removed.
  The library will now be compiled with the default C++ standard enabled by
  the compiler. This is (as of May 2020) C++14 for all compilers. If you want
  to override that behavior, please set the C++ standard directly for example
  by configuring with <code>-DDEAL_II_CXX_FLAGS="-std=c++17"</code>, or by
  setting the environment variable <code>CXXFLAGS="-std=c++17"</code>.
  <br>
  (Matthias Maier, 2020/05/21)
 </li>

 <li>
  Updated: deal.II now requires a compiler with enabled C++14 support.
  <br>
  (Matthias Maier, 2020/05/21)
 </li>

 <li>
  Changed: The polynomial space template argument has been removed from
  FE_Poly and FE_PolyTensor.
  <br>
  (Graham Harper, Daniel Arndt, 2020/05/21)
 </li>

 <li>
  Removed: The deprecated class Threads::PosixThreadBarrier has been
  removed.
  <br>
  (Wolfgang Bangerth, 2020/04/21)
 </li>

</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="920-930-general"></a>
<h3>General</h3>
<ol>

 <li>
  New: The step-78 tutorial program solves the Black-Scholes equations.
  <br>
  (Tylor Anderson, 2021/05/19)
 </li>

 <li>
  New: The step-79 tutorial program demonstrates solution of a topology optimization problem.
  <br>
  (Justin O'Connor, 2021/05/19)
 </li>

 <li>
  New: The step-66 tutorial program shows how to solve a nonlinear PDE in parallel
  with MatrixFree methods.
  <br>
  (Fabian Castelli, 2021/05/19)
 </li>

 <li>
  New: The class FEPointEvaluation provides an interface to the evaluation of
  interpolated solution values and gradients on cells on arbitrary reference
  point positions. These points can change from cell to cell, both with
  respect to their quantity as well to the location. The two typical use
  cases are evaluations on non-matching grids and particle simulations.
  <br>
  (Martin Kronbichler, 2021/05/18)
 </li>

 <li>
  New: The step-72 tutorial revisits step-15 (that solves a nonlinear elliptic
  problem using Newton's method). In this variant, automatic differentiation
  is used either to linearize the residual, or to compute the full linear system
  from an energy functional.
  <br>
  (Jean-Paul Pelteret, 2021/05/06)
 </li>

 <li>
  New: The step-71 tutorial demonstrates how automatic and symbolic differentiation
  can be leveraged to compute the derivatives of a coupled constitutive law.
  <br>
  (Jean-Paul Pelteret, 2021/05/06)
 </li>

 <li>
  New: The step-77 tutorial explores one of the mentioned extensions to step-15,
  namely implementing a proper step length control algorithm. It uses the KINSOL
  nonlinear solver package from the SUNDIALS library to incorporate a line
  search method into the incremental update for Newton's method.
  <br>
  (Wolfgang Bangerth, 2021/05/04)
 </li>

 <li>
  New: The step-75 tutorial program shows how to solve a simple Laplace
  equation in parallel with hp-adaptive and MatrixFree methods.
  <br>
  (Marc Fehling, Peter Munch, Wolfgang Bangerth, 2021/04/15)
 </li>

 <li>
  New: Experimental support for simplex meshes (i.e., meshes consisting
  of triangle or tetrahedron cells) and mixed meshes (i.e., meshes consisting
  of lines, triangles, quadrilateral, tetrahedron, pyramid, wedge, and/or
  hexahedron cells) has been added. The new features are presented
  in the test folder `tests/simplex` and in the
  new "Simplex" module page.
  <br>
  (various, 2021/03/04)
 </li>

 <li>
  New: The class Utilities::MPI::RemotePointEvaluation and the function
  VectorTools::point_values() allow to work on arbitrary distributed
  points.
  <br>
  (Peter Munch, Martin Kronbichler, Magdalena Schreter, Niklas Fehn, 2021/02/28)
 </li>

 <li>
  New: The communicator of an arbitrary (not just parallel) Triangulation class can now be
  queried via Triangulation::get_communicator() or DoFHandler::get_communicator(). In
  the case of serial Triangulations and DoFHandler set up with serial Triangulations,
  MPI_COMM_SELF is returned.
  <br>
  (Peter Munch, 2021/02/28)
 </li>

 <li>
  New: The behavior of the local_size() member function is not consistent across
  all vector classes that support ghost elements. As a remedy this member
  function is deprecated and replaced by locally_owned_size() that returns the
  number of locally owned elements (in particular without ghost elements).
  <br>
  (Daniel Arndt, David Wells, 2021/02/11)
 </li>

 <li>
  New: Added a new quadrature rule QWitherdenVincent for simplices.
  <br>
  (David Wells, 2021/02/08)
 </li>

 <li>
  New: Added a new class BarycentricPolynomial that makes defining
  polynomials on simplices much easier.
  <br>
  (David Wells, 2021/01/26)
 </li>

 <li>
  New: Added a new finite element FE_SimplexP_Bubbles suitable for using mass
  lumping on simplex meshes.
  <br>
  (David Wells, 2021/01/26)
 </li>

 <li>
  New: The step-74 tutorial implements the symmetric interior penalty Galerkin
  (SIPG) method for Poisson's equation using the FEInterfaceValues class (for
  interface terms) in conjunction with the MeshWorker::mesh_loop() concept.
  <br>
  (Timo Heister, Jiaqi Zhang, 2021/01/04)
 </li>

 <li>
  Improved: Update SUNDIALS ARKODE interface to support versions > 4.0.0 and add preconditioner support
  <br>
  (Sebastian Proell, 2020/12/11)
 </li>

 <li>
  New: Added capacity to update ghost particles without rebuilding them from scratch in the particle_handler
  <br>
  (Bruno Blais, Peter Munch, 2020/11/10)
 </li>

 <li>
  New: Tutorial example (step-68) showcasing parallel simulation of the advection
  of particles including load balancing.
  <br>
  (Bruno Blais, Toni El Geitani Nehme, Rene Gassm
  ller, Peter Munch, 2020/05/23)
 </li>

 <li>
  New: The step-19 tutorial program shows how to use particle methods.
  <br>
  (Wolfgang Bangerth, Rene Gassmoeller, Peter Munch, 2020/09/15)
 </li>

 <li>
  New: The GridIn class can now read ExodusII files when deal.II is configured
  with Trilinos and SEACAS.
  <br>
  (David Wells, 2020/09/09)
 </li>

 <li>
  Changed: The internal data structures of DoFHandler have been modified to use
  compressed row storage, enabling it to also handle hp::DoFHandler functionalities.
  Currently, the user can choose between the normal mode and the hp mode during
  calling the constructor. Please note that the multigrid functionalities are only
  available during normal mode.
  <br>
  (Peter Munch, 2020/05/23)
 </li>

</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="920-930-specific"></a>
<h3>Specific improvements</h3>
<ol>

 <li>
  Fixed: Make SUNDIALS::IDA function names  consistent with SUNDIALS::KINSOL function names.
  <br>
  (Luca Heltai, 2021/05/21)
 </li>

 <li>
  Fixed: KINSOL interface is now fully working
  <br>
  (Luca Heltai, Nicola Giuliani, 2021/05/19)
 </li>

 <li>
  Improved: The images in tutorial 35 have been updated to cold and hot color scheme with higher resolution and more consistent sizing.
  <br>
  (Bryn Barker, 2021/05/18)
 </li>

 <li>
  New: SUNDIALS::IDA is now compatible with SUNDIALS > 4.0.0
  <br>
  (Nicola Giuliani, Luca Heltai, 2021/05/13)
 </li>

 <li>
  Improved: ParticleHandler::sort_particles_into_subdomains_and_cells()
  now uses a faster method to identify particles in cells farther away
  from their previous cell.
  <br>
  (Martin Kronbichler, 2021/05/12)
 </li>

 <li>
  Fixed: DataOutBase::write_ucd() has been fixed to output vertices in the correct
  order. Previously, cells were tangled when read with Paraview's AVS-UCD reader.
  <br>
  (David Wells, 2021/05/06)
 </li>

 <li>
  Fixed: The function Particles::Generators::probabilistic_locations could crash
  if some MPI ranks had no active cells. This is fixed now.
  <br>
  (Rene Gassmoeller, 2021/05/05)
 </li>

 <li>
  Improved: There is now a cache for TriaAccessor::vertex_index() for the cells
  of a mesh, which considerably speeds up operations with intensive use of
  `cell->vertex()` or `cell->vertex_index()`.
  <br>
  (Martin Kronbichler, 2021/05/04)
 </li>

 <li>
  New: Add a function taylor_estimate_function_bounds, which estimates the range
  of the value and gradient components of a Function<dim> over a
  BoundingBox<dim>, by approximating the function by a 2nd order Taylor
  polynomial.
  <br>
  (Simon Sticko, 2021/04/25)
 </li>

 <li>
  Deprecated: The FEValuesViews::Scalar, FEValuesViews::Vector,
  FEValuesViews::Tensor, and FEValuesViews::SymmetricTensor classes all
  have `OutputType` structures that defined the types of evaluating
  the values and derivatives of finite element fields at quadrature
  points. These structures have been deprecated, and the corresponding
  types are now defined via local `using` declarations in the classes
  mentioned above.
  <br>
  (Wolfgang Bangerth, 2021/04/25)
 </li>

 <li>
  Bugfix: FE_Nedelec<2>::convert_generalized_support_point_values_to_dof_values()
  now works correctly for every degree.
  <br>
  (Jake Harmon, 2021/04/22)
 </li>

 <li>
  Fixed: The function hp::Refinement::predict_error() produced incorrect
  results for p-coarsening.
  <br>
  (Marc Fehling, 2021/04/21)
 </li>

 <li>
  Deprecated: The TriaAccessor::number_of_children() function has been
  deprecated in favor of the new TriaAccessor::n_active_descendants()
  function.
  <br>
  (Wolfgang Bangerth, 2021/04/20)
 </li>

 <li>
  Added: The MPI_InitFinalize RAII class has gained an
  MPI_InitFinalize::signals::at_mpi_init and an
  MPI_InitFinalize::signals::at_mpi_finalize signal that are triggered
  immediately after initializing the MPI context with <code>MPI_Init</code>
  and immediately before deinitializing the MPI context with
  <code>MPI_Finalize</code>.
  <br>
  (Matthias Maier, 2021/04/19)
 </li>

 <li>
  New: There are now versions of the constructors of the
  Functions::InterpolatedTensorProductGridData and
  Functions::InterpolatedUniformGridData classes
  that *move* the data in the arguments given,
  instead of copying it.
  <br>
  (Wolfgang Bangerth, 2021/04/19)
 </li>

 <li>
  Fixed: Work around a memory leak issue in OpenMPI 4.1.0 triggered by our
  Utilities::MPI::min_max_avg() function by repeatedly allocating and freeing
  MPI_Datatype handles.
  <br>
  (Matthias Maier, 2021/04/18)
 </li>

 <li>
  New: FEInterfaceValues now can use FEValuesExtractors to extract scalar or
  vector components like FEValues does.
  <br>
  (Jiaqi Zhang, 2021/04/15)
 </li>

 <li>
  New: The functions AlignedVector::replicate_across_communicator()
  and  Table::replicate_across_communicator()
  allow copying information across MPI processes.
  <br>
  (Wolfgang Bangerth, 2021/04/28)
 </li>

 <li>
  New: Now the HDF5 interface can set bool attributes.
  <br>
  (Daniel Garcia-Sanchez, 2021/04/03)
 </li>

 <li>
  Add a new grid generator (subdivided_cylinder) which generates a
  cylinder with a number of x subdivisions which is specified by the user.
  Alter cylinder so that it also uses this new generator to prevent code
  duplication.
  <br>
  (Bruno Blais, 2021/04/01)
 </li>

 <li>
  Add: Adds new settings to PETScWrappers::PreconditionBoomerAMG::initialize()
  in the struct PETScWrappers::PreconditionBoomerAMG::AdditionalData.
  <br>
  (Maximilian Bergbauer, 2021/03/31)
 </li>

 <li>
  Deprecated: The version of DoFTools::extract_boundary_dofs() that
  returns its information via an `IndexSet` reference argument has been
  deprecated. Use the version of the function that returns information
  via an IndexSet return type instead.
  <br>
  (Wolfgang Bangerth, 2021/03/30)
 </li>

 <li>
  Deprecated: The version of DoFTools::extract_boundary_dofs() that
  returns its information via a `std::vector<bool>` has been
  deprecated. Use the version of the function that returns information
  via an IndexSet instead.
  <br>
  (Wolfgang Bangerth, 2021/03/30)
 </li>

 <li>
  New: Class parallel::distributed::TemporarilyMatchRefineFlags that
  temporarily modifies the refine and coarsen flags of all active cells
  on a parallel::distributed::Triangulation to match its p4est oracle.
  <br>
  (Marc Fehling, 2021/03/29)
 </li>

 <li>
  New: The getter function for the divergence in
  FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_divergence()
  is now also implemented for dim==1 and n_components==1. The function
  FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_gradient()
  now also takes a rank 2 tensor as input argument.
  <br>
  (Magdalena Schreter, Peter Munch, Martin Kronbichler, 2021/03/22)
 </li>

 <li>
  New: Add nearest neighbor predicate support to ArborXWrappers::BVH(). This
  allows to find the `n` closest bounding boxes or points to any given bounding
  box or point.
  <br>
  (Bruno Turcksin, 2021/03/19)
 </li>

 <li>
  Improved: The function parallel::distributed::Triangulation::copy_triangulation()
  now also works with locally refined meshes of type
  parallel::distributed::Triangulation.
  <br>
  (Marc Fehling, 2021/03/19)
 </li>

 <li>
  Improved: The class LinearAlgebra::ReadWriteVector now also can import
  from Vector and LinearAlgebra::Vector.
  <br>
  (Peter Munch, 2021/03/18)
 </li>

 <li>
  New: Adds the new generalized MPI_Bcast function Utilities::MPI::broadcast() for arbitrary data types T.
  <br>
  (Maximilian Bergbauer, 2021/03/18)
 </li>

 <li>
  Fixed: MatrixFree::get_boundary_id() would previously only support boundary ids
  between 0 and 255 via an `unsigned char` type. This has now been changed to
  the correct types::boundary_id type, supporting also larger boundary ids.
  <br>
  (Martin Kronbichler, 2021/03/16)
 </li>

 <li>
  New: Added Functions::IdentityFunction.
  <br>
  (David Wells, 2021/03/16)
 </li>

 <li>
  New: The new method Differentiation::SD::BatchOptimizer::extract() allows one to
  extract results from previous evaluation of the symbolic expressions.
  This means that one can safely use a single instance of a batch optimizer to
  pre-compute a set of results, cache them, and then later extract some results
  from the cached data without re-evaluating any symbolic expressions.
  <br>
  (Jean-Paul Pelteret, 2021/03/13)
 </li>

 <li>
  Improved: The result type deduction for the product of symbolic types
  (specifically, Differentiation::SD::Expression) with tensors and symmetric
  tensors has been improved.
  <br>
  (Jean-Paul Pelteret, 2021/03/13)
 </li>

 <li>
  New: The MeshWorker::ScratchData class is now able to compute the Laplacians
  of the scalar-valued and vector-valued components of the solution vector.
  <br>
  (Jean-Paul Pelteret, 2021/03/13)
 </li>

 <li>
  Improved: The MeshWorker::ScratchData::get_general_data_storage() method now has
  a `const` variant.
  <br>
  (Jean-Paul Pelteret, 2021/03/13)
 </li>

 <li>
  Improved: The function VectorTools::get_position_vector() can now also take
  a Mapping object as input parameter.
  <br>
  (Peter Munch, 2021/03/10)
 </li>

 <li>
  New: MappingQCache has new initialize functions. One set of functions takes
  either a std::function or a dealii::Function for transforming
  individual points. Another set of functions take global vectors (and
  MGLevelObjects of global vectors) and use these to initialize the position of
  the support points much like MappingFEField does.
  <br>
  (Niklas Fehn, Martin Kronbichler, Peter Munch, 2021/03/10)
 </li>

 <li>
  New: Utilities::MPI::logical_or() for collective <i>logical or</i> operations.
  <br>
  (Marc Fehling, 2021/03/04)
 </li>

 <li>
  New: Add a component parameter to FESeries::Fourier/Legendre to
  make them working with non-primitive elements.
  <br>
  (Ce Qin, 2021/03/01)
 </li>

 <li>
  New: A new class GridTools::MarchingCubeAlgorithm has been added. It
  helps to create a surface mesh on the iso line/contour of a scalar field.
  <br>
  (Peter Munch, Magdalena Schreter, Martin Kronbichler, 2021/02/28)
 </li>

 <li>
  New: Function hp::Refinement::limit_p_level_difference()
  restricts the maximal level difference of neighboring cells with respect
  to the finite element hierarchy of the registered hp::FECollection.
  <br>
  (Marc Fehling, 2021/02/23)
 </li>

 <li>
  New: The new MGTransferBase::prolongate_and_add() performs a prolongation
  without zeroing out the destination vector.
  <br>
  (Martin Kronbichler, 2021/02/22)
 </li>

 <li>
  New: Member function hp::FECollection::get_hierarchy_sequence() returning
  the sequence of finite element indices that correspond to the registered
  hierarchy.
  <br>
  (Marc Fehling, 2021/02/22)
 </li>

 <li>
  Fixed: deal.II should now be able to use versions of LAPACK compiled with
  Fortran compilers that do not adhere to the usual Fortran mangling
  convention.
  <br>
  (David Wells, 2021/02/20)
 </li>

 <li>
  Added a template version of CUDAWrappers::MatrixFree::reinit() that takes an
  IteratorFilters object. This allows to perform the operator evaluation on part
  of the domain.
  <br>
  (Bruno Turcksin, 2021/02/19)
 </li>

 <li>
  Added: A multidimensional array dealii::ndarray that allows one to
  conveniently create "stacked" std::array objects that model a
  multidimensional array. `dealii::ndarray<double, 1, 2, 3, 4>` for example
  is a short-hand for `std::array<std::array<std::array<std::array<double, 4>, 3>, 2>, 1>`.
  <br>
  (Matthias Maier, 2021/02/19)
 </li>

 <li>
  Moved: The CommunicationPatternBase class has been moved from the LinearAlgebra
  namespace into the more general Utilities::MPI namespace.
  <br>
  (David Wells, 2021/02/15)
 </li>

 <li>
  Fixed: CutOffFunctionTensorProduct::gradient() was computing the wrong gradient. This is now fixed.
  <br>
  (Luca Heltai, 2021/02/13)
 </li>

 <li>
  New: ArrayView objects can now also be constructed from C-style arrays.
  <br>
  (Wolfgang Bangerth, 2021/02/09)
 </li>

 <li>
  New: The method parallel::distributed::Triangulation::load()
  can now also accept a p4est/p8est forest, which can be queried
  from an existing triangulation via parallel::distributed::Triangulation::get_p4est().
  <br>
  (Marc Fehling, Peter Munch, 2021/02/09)
 </li>

 <li>
  New: There are now functions Utilities::get_bit() and
  Utilities::set_bit() that do as their names suggest.
  <br>
  (Peter Munch, Wolfgang Bangerth, 2021/02/05)
 </li>

 <li>
  Improved: AffineConstraints::copy_from() now also works for differing
  number template types.
  <br>
  (Peter Munch, Maximilian Bergbauer, 2021/02/05)
 </li>

 <li>
  Update: class PETScWrappers::PreconditionerBase renamed to
  PETScWrappers::PreconditionBase for consistency.
  <br>
  (Pasquale Claudio Africa, 2021/02/03)
 </li>

 <li>
  Fixed: Triangulation::get_manifold_ids() was not returning all ids correctly.
  <br>
  (Luca Heltai, 2021/02/03)
 </li>

 <li>
  New: SUNDIALS N_Vector module now directly operates on different deal.II vectors without internally
  creating copies.
  In particular, ARKode can now also be used with LinearAlgebra::distributed::(Block)Vector.
  <br>
  (Sebastian Proell, 2021/01/31)
 </li>

 <li>
  Improved: Exception texts are now formatted and broken to fixed-length
  lines, rather than flowing to arbitrary lengths.
  <br>
  (Wolfgang Bangerth, 2021/01/25)
 </li>

 <li>
  New: Added support for gmsh library API. This allows using GridIn::read_msh()
  and GridOut::write_msh() to save and read also manifold id information.
  <br>
  (Luca Heltai, 2021/01/20)
 </li>

 <li>
  New: implemented p::f::Triangulation::load()/save() for use with p::d::SolutionTransfer.
  <br>
  (Pasquale Claudio Africa, Peter Munch, 2021/01/18)
 </li>

 <li>
  Improved: provide abstract interface to p::d::Triangulation::load()/save().
  <br>
  (Pasquale Claudio Africa, 2021/01/18)
 </li>

 <li>
  New: Added python wrappers to enable support for simplex and mixed meshes.
  <br>
  (Alexander Grayver, 2021/01/18)
 </li>

 <li>
  New: Mapping::get_bounding_box() now returns the bounding box of the support points, when the mapping
  is of type MappingQ, MappingQGeneric, and MappingQEulerian.
  <br>
  (Luca Heltai, 2021/01/16)
 </li>

 <li>
  New: The method parallel::distributed::Triagnulation::communicate_locally_moved_vertices()
  has been refactored and moved to parallel::DistributedTriangulationBase so that it can now also be
  used for parallel::fullydistributed::Triangulation.
  <br>
  (Daniel Arndt, Ivan Fumagalli, Peter Munch, 2021/01/15)
 </li>

 <li>
  New: Implemented MappingFE::transform_real_to_unit_cell().
  <br>
  (Luca Heltai, 2021/01/13)
 </li>

 <li>
  New: The old tensor basis transformation functions internal::Physics::transformation_contraction()
  have been moved out of the internal namespace and renamed to
  Physics::Transformations::basis_transformation() and have documentation now.
  <br>
  (Nils Much, 2021/01/12)
 </li>

 <li>
  New: Created two new methods ParameterAcceptor::enter_subsection() and
  ParameterAcceptor::leave_subsection() that allow more intricated parameter
  paths within ParameterAcceptor classes.
  <br>
  (Luca Heltai, 2021/01/12)
 </li>

 <li>
  New: Add new BVH (Bounding Volume Hierarchy) class based on ArborX. This class
  performs multiple kinds of geometric search: intersection of bounding boxes and
  intersection of bounding boxes with points.
  <br>
  (Bruno Turcksin, 2021/01/08)
 </li>

 <li>
  New: Function GridTools::get_subdomain_association() determines
  the owning process of any valid cell on a Triangulation that is
  represented by a globally unique CellId.
  <br>
  (Marc Fehling, 2021/01/07)
 </li>

 <li>
  New: The policy under which things in deal.II are deprecated has changed.
  Deprecated features are now first marked with DEAL_II_DEPRECATED_EARLY until the
  next release of the library, at which point they will be remarked with
  DEAL_II_DEPRECATED. By default, things marked with DEAL_II_DEPRECATED_EARLY do
  not print deprecation warnings - this is controlled with the
  DEAL_II_EARLY_DEPRECATIONS CMake configuration option. This change was made so
  that users can use multiple recent checkouts of the development branch without
  needing to address the problem that some will print deprecation warnings and
  others do not, and also so that new deprecation warnings do not appear outside
  of the release period.
  <br>
  (Daniel Arndt, 2021/01/05)
 </li>

 <li>
  Fixed: Reset time and timestep_number during pre refinement steps in step-26.
  <br>
  (Praveen Chandrashekar, 2021/01/05)
 </li>

 <li>
  Changed: Tutorial step-27 has been simplified and now uses the recently
  introduced SmoothnessEstimator namespace.
  <br>
  (Marc Fehling, 2020/12/24)
 </li>

 <li>
  Fixed: DoFTools::make_periodicity_constraints() was not instantiated for condition 'dim < spacedim'. The instantiation was corrected and a test was added to verify that the function works for such condition.
  <br>
  (Malhar Tidke, 2020/12/23)
 </li>

 <li>
  Improved: If there is an uncaught exception the destructor of the HDF5
  interface does not call H5Dclose, H5Gclose or H5Fclose. In addition,
  to avoid MPI synchronization and a possible deadlock, the destructor
  calls MPI_Abort().
  (Daniel Garcia-Sanchez, 2020/12/23)
 </li>

 <li>
  New: Add functions CUDAWrappers::local_q_point_id_host(),
  CUDAWrappers::get_quadrature_point_host(), CUDAWrappers::copy_mf_data_to_host(), and
  CUDAWrappers::MatrixFree::get_colored_graph(). These functions can be used to
  evaluate material properties in the same order on the host and on the device.
  <br>
  (Bruno Turcksin, 2020/12/06)
 </li>

 <li>
  Fixed: A previous patch accidentally broke step-46 and led to
  exceptions about accessing neighboring cells that don't actually exist
  because a cell is at the boundary. This has been fixed.
  <br>
  (Wolfgang Bangerth, 2020/12/03)
 </li>

 <li>
  New: The function GridGenerator::convert_hypercube_to_simplex_mesh allows to
  convert a given triangulation based on quadrilaterals (2D) or hexahedra (3D)
  to a triangulation based on simplices or tetraedra, respectively. Thereby,
  material_IDs and boundary_IDs are inherited from the given triangulation.
  <br>
  (Elias Dejene, Peter Munch, 2020/11/23)
 </li>

 <li>
  New: A new class LowStorageRungeKutta is added to the namespace TimeStepping to
  implement the explicit low-storage Runge-Kutta methods, see @cite KennedyCarpenterLewis2000 and step-67.
  <br>
  (Jiaqi Zhang, 2020/11/18)
 </li>

 <li>
  Improved: QCollection<dim> now also accepts Quadrature<1> and converts this
  input quadrature rule to a dim-dimensional quadrature rule internally via
  tensor product. Furthermore, a copy constructor has been added
  accepting QCollection<1>.
  <br>
  (Peter Munch, 2020/11/12)
 </li>

 <li>
  Fixed: WorkStream::mesh_loop() should now work on anisotropic grids.
  <br>
  (Luca Heltai, 2020/11/03)
 </li>

 <li>
  Improved: An inverse quadratic approximation has been added for the pull-back
  operation in the TransfiniteInterpolationManifold::new_points() function. The
  better initial guesses for the Newton/Broyden iteration make the computation
  faster and more robust in some difficult scenarios.
  <br>
  (Martin Kronbichler, 2020/11/01)
 </li>

 <li>
  Fixed: There was a bug in
  Differentiation::AD::ScalarFunction::extract_hessian_component()
  that was triggered when using a symmetric tensor extractor with a non-zero
  first component index. Having previously lead to either incorrect results being
  returned or valid user programs crashing, it has now been corrected.
  <br>
  (Jean-Paul Pelteret, 2020/11/01)
 </li>

 <li>
  Improved: MatrixFree now also works for hp in MPI-parallelized
  programs.
  <br>
  (Marc Fehling, Katharina Kormann, Martin Kronbichler, Peter Munch, 2020/10/20)
 </li>

 <li>
  New: step-9 uses the "streamline-upwind Petrov-Galerkin" method, but
  does not make any attempt at explaining what this method is or why it
  might be named like this. This has been rectified: The introduction
  now has a long section that explains the origin of the method and its name.
  <br>
  (Wolfgang Bangerth, 2020/10/10)
 </li>

 <li>
  New: Mapping::transform_points_real_to_unit_cell() can compute the operation
  of Mapping::transform_real_to_unit_cell() on many points simultaneously, which
  can be much faster for MappingQGeneric and derived classes that involve
  expensive operations to compute the support points of the mapping.
  <br>
  (Martin Kronbichler, 2020/10/07)
 </li>

 <li>
  New: GridRefinement::refine_and_coarsen_fixed_fraction() and
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction()
  now allow to specify a VectorTools::NormType, which determines how
  combined errors on subsets of cells will be calculated.
  <br>
  (Marc Fehling, 2020/10/06)
 </li>

 <li>
  New: SSP_THIRD_ORDER is added to the namespace TimeStepping to
  implement the explicit third order Strong Stability Preserving (SSP) Runge-Kutta method,
  which is also called the third order Total Variation Diminishing (TVD) Runge-Kutta method, see @cite gottlieb2001strong.
  <br>
  (Jiaqi Zhang, 2020/10/05)
 </li>

 <li>
  New: GridTools::affine_cell_approximation() returns a matrix <i>A</i> and
  offset vector <i>b</i> that describe a least-squares fit of an affine
  approximation to a set of vertices of a cell.
  <br>
  (Martin Kronbichler, 2020/10/04)
 </li>

 <li>
  New: Helper functions CellAccessor::child_iterators() and
  DoFCellAccessor::child_iterators() which return iterators to children of
  a cell via `cell->child_iterators()`.
  <br>
  (Marc Fehling, 2020/10/03)
 </li>

 <li>
  New: CellId has a new constructor to create it from a std::string.
  <br>
  (Timo Heister, 2020/10/05)
 </li>

 <li>
  Improved: MappingQGeneric::transform_real_to_unit_cell() has been made much
  faster by directly working with the tensor product form of the mapping shape
  functions and avoiding many unnecessary memory allocations. The main cost is
  now MappingQGeneric::compute_mapping_support_points(), which can be made fast
  with MappingQCache, for example.
  <br>
  (Martin Kronbichler, 2020/09/30)
 </li>

 <li>
  New: The function BlockSparsityPattern::print_svg() outputs a block
  sparsity pattern in SVG format.
  <br>
  (Wolfgang Bangerth, 2020/09/25)
 </li>

 <li>
  Changed: step-29 no longer uses the `deallog` variable to generate
  output, but instead directly writes to `std::cout`.
  <br>
  (Wolfgang Bangerth, 2020/09/23)
 </li>

 <li>
  New: The classes FEEvaluation and FEFaceEvaluation with template parameter -1
  for the polynomial degree is now based on pre-compiled templated code for
  polynomial degrees between 1 and 6. This allows for fast execution of
  matrix-free evaluation for run-time polynomial degrees. The generated
  instantiations are controlled by
  `include/deal.II/matrix_free/evaluation_template_factory.templates.h` and can
  be pre-compiled for additional degrees in user code.
  <br>
  (Martin Kronbichler, Peter Munch, 2020/09/21)
 </li>

 <li>
  Fixed: Our cmake scripts forgot to install some of the files that are
  part of the code gallery. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/09/17)
 </li>

 <li>
  Fixed: DoFTools::extract_dofs() return an IndexSet as result used to have a
  quadratic complexity in the number of extracted indices. This is now fixed.
  <br>
  (Martin Kronbichler, 2020/08/11)
 </li>

 <li>
  Added: method for returning list of all triangulation cells.
  <br>
  (Alexander Grayver, 2020/09/03)
 </li>

 <li>
  Added: python wrapper for GridTools::replicate_triangulation,
  more general version of the GridTools::merge_triangulations is
  implemented
  <br>
  (Alexander Grayver, 2020/09/01)
 </li>

 <li>
  New: The methods FEEvaluation::gather_evaluate(),
  FEEFaceEvaluation::gather_evaluate(), FEEvaluation::integrate_scatter() and
  FEfaceEvaluation::integrate_scatter() can now also accept block vectors.
  <br>
  (Peter Munch, Magdalena Schreter, Martin Kronbichler, 2020/08/31)
 </li>

 <li>
  New: A particle collection can now be copied into a new ParticleHandler object using
  the new ParticleHandler::copy_from function.
  <br>
  (Rene Gassmoeller, 2020/08/28)
 </li>

 <li>
  Added: DataOut now supports HDF5 file format with simplex meshes.
  <br>
  (Pasquale Claudio Africa, 2020/08/27)
 </li>

 <li>
  Updated: GridTools::transform() now works with simplex meshes.
  <br>
  (Pasquale Claudio Africa, 2020/08/26)
 </li>

 <li>
  Improved: The macro DEAL_II_PICKUP_TESTS can now also be run with on optional
  input parameter that can be used to manipulate the folder name of tests during
  ctest.
  <br>
  (Peter Munch, 2020/07/23)
 </li>

 <li>
  New: The test suite can now also be run with .json files.
  <br>
  (Peter Munch, 2020/08/17)
 </li>

 <li>
  Improved: The definitions of the templated functions of the HDF5 interface are
  now in hdf5.h, therefore it is possible to use additional template combinations.
  The instantiations are no longer necessary, therefore they have been removed.
  <br>
  (Daniel Garcia-Sanchez, 2020/08/14)
 </li>

 <li>
  New: MGTools::make_sparsity_pattern() can now take an optional
  AffineConstraints argument to add the effect of, e.g., periodic boundary
  conditions.
  <br>
  (Martin Kronbichler, 2020/08/11)
 </li>

 <li>
  Fixed: The DataPostprocessorTensor class erroneously announced that
  the components of its output are a bunch of scalars when, of course,
  the whole point of the class was to output things as one tensor
  field. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/08/06)
 </li>

 <li>
  New: GridIn::read_vtk() now supports triangular and tetrahedral meshes.
  <br>
  (Peter Munch, 2020/07/23)
 </li>

 <li>
  New: GridIn::read_msh() now supports triangular and tetrahedral meshes.
  <br>
  (Daniel Paukner, 2020/07/20)
 </li>

 <li>
  New: The method hp::FEFaceValues::reinit() can now also accept face iterators.
  <br>
  (Peter Munch, 2020/07/16)
 </li>

 <li>
  New: The class hp::MappingCollection has a new constructor. This constructor creates
  a MappingCollection from one or more mapping objects passed to the constructor.
  <br>
  (Peter Munch, 2020/07/15)
 </li>

 <li>
  Fixed: MeshWorker::mesh_loop() did not work on 1d with refined grids. This is now fixed.
  <br>
  (Luca Heltai, 2020/07/08)
 </li>

 <li>
  Added the functions CUDAWrappers::MatrixFree::get_vector_partitioner() and
  CUDAWrappers::MatrixFree::get_dof_handler()
  <br>
  (Bruno Turcksin, 2020/07/06)
 </li>

 <li>
  Fixed: In parallel hp-adaptive applications,
  DoFHandler::distribute_dofs() no longer fails to enumerate degrees of
  freedom on ghost interfaces if continuous finite elements do not
  dominate each other.
  <br>
  (Marc Fehling, 2020/07/03)
 </li>

 <li>
  New: A new quadrature rule for simplex geometric entities has been added.
  <br>
  (Peter Munch, 2020/07/02)
 </li>

 <li>
  New: Geometric objects of a Triangulation are assigned a ReferenceCell::Type. The value
  can be queried via TriaAccessor::reference_cell_type().
  <br>
  (Peter Munch, 2020/06/30)
 </li>

 <li>
  New: The class ArrayView can now also be constructed from std::array.
  <br>
  (Peter Munch, 2020/06/29)
 </li>

 <li>
  New: BoundingBox::real_to_unit() and BoundingBox::unit_to_real() allow one to
  apply both the direct and the inverse transformation that are needed to map the
  unit bounding box to the current box, and vice-versa.
  <br>
  (Luca Heltai, 2020/06/29)
 </li>

 <li>
  New: The member function DiscreteTime::set_next_step_size() is added.
  <br>
  (Reza Rastak, 2020/06/27)
 </li>

 <li>
  New: There is now a constructor for class Tensor that takes
  an initializer from an object of type ArrayView.
  <br>
  (Wolfgang Bangerth, 2020/06/27)
 </li>

 <li>
  Fixed: The class parallel::distributed::SolutionTransfer can now
  also handle FE_Nothing.
  <br>
  (Dominic Soldner, Peter Munch, 2020/06/24)
 </li>

 <li>
  Fixed: FEInterfaceValues now works also for codim one and two. Instantiated
  also DoFTools::make_flux_sparsity_pattern() for codim one and two.
  <br>
  (Luca Heltai, 2020/06/24)
 </li>

 <li>
  New: The class TriaAccessor provides now the capability to query
  the number of vertices, lines, and faces (with n_vertices(),
  n_lines(), n_faces(), vertex_indices(),
  line_indices(), face_indices()). The new methods can be used as an
  alternative to the methods in GeometryInfo.
  <br>
  (Peter Munch, 2020/06/23)
 </li>

 <li>
  New: Added FEInterfaceValues to MeshWorker::ScratchData.
  <br>
  (Luca Heltai, 2020/06/23)
 </li>

 <li>
  Fixed: The ParticleHandler::insert_particles() function forgot to
  associate particles with the common property pool. Consequently,
  setting properties on particles added to a ParticleHandler this way
  led to an error.
  <br>
  (Andrew Davis, Wolfgang Bangerth, 2020/06/23)
 </li>

 <li>
  New: Particles::Particle and Particles::ParticleAccessor can now be used as
  indexable in boost::rtree objects.
  <br> (Luca Heltai, 2020/06/15)
 </li>

 <li>
  New: Each cell is assigned a globally unique active cell index and (if requested)
  a level cell index. These indices are integers enumerated contiguously within
  each subdomain of the mesh.
  Users can query locally-owned and ghost cells for their indices via CellAccessor::global_active_cell_index()
  or CellAccessor::global_level_cell_index().
  The value is managed automatically by the Triangulation classes.
  Furthermore, triangulations deriving from parallel::TriangulationBase provide partitioners
  for these indices, which can be used to set up ghosted vectors with one entry per
  cell.
  <br>
  (Peter Munch, 2020/06/12)
 </li>

 <li>
  New: The function Particles::ParticleHandler::add_global_particles() now takes
  another optional argument, that allows one to set ids arbitrarily. Moreover,
  now the numbering of the ids is correct also if we call the method more than
  one time. Newly added particles, if ids are not specified, now correctly get
  the first available ids.
  Added a new version of Particles::ParticleHandler::add_global_particles() that
  takes a vector of Particles::Particle objects instead of just their positions.
  This can be used in conjunction with the signal
  Particles::ParticleHandler::Signals::particle_lost() to reinsert
  Particles::Particle objects that went out of the locally owned and ghost cells.
  <br> (Luca Heltai, 2020/06/11)
 </li>

 <li>
  Fixed: The AffineConstraints class had a bug where, if deal.II was
  compiled without threading, the class used a global variable. If a
  user program used threads nonetheless, this global variable led to
  race conditions. This is now fixed.
  <br>
  (Wolfgang Bangerth, 2020/06/11)
 </li>

 <li>
  Fixed: Fix a bug where CUDAWrappers::MatrixFree::cell_loop() would set the
  destination vector to zero if the partitioner of the MatrixFree object was
  different from the partitioner of the source or destination vector.
  <br>
  (Bruno Turcksin, 2020/06/10)
 </li>

 <li>
  New: Applying user constraints
  before prolongation in MGTransferPrebuilt.
  <br>
  (Julian Roth and Conrad Clevenger, 2020/06/05)
 </li>

 <li>
  New: FEEvaluation::evaluate(), FEEvaluation::integrate(),
  FEFaceEvaluation::evaluate() and FEFaceEvaluation::integrate() now take an
  argument of type EvaluationFlags::EvaluationFlags to determine which of the
  values, gradients or hessians should be evaluated to quadrature points or
  integrated, allowing much more expressive programs than the previous list of
  bools. The evaluation flags can be combined with `operator|`, similarly to
  UpdateFlags for FEValues.
  <br>
  (Timo Heister, 2020/06/05)
 </li>

 <li>
  Changed: The vertices in CellData are now stored in form of a std::vector
  instead of C-style array.
  <br>
  (Peter Munch, 2020/05/31)
 </li>

 <li>
  Improved: The efficiency of the assembly of step-62 has been improved and now
  it is 7 times faster.
  <br>
  (Daniel Garcia-Sanchez, 2020/05/31)
 </li>

 <li>
  Fixed: Fix a bug where only one CUDAWrappers::MatrixFree object was valid at a
  given time. There is now a variable CUDAWrappers::mf_n_concurrent_objects in
  base/cuda_size.h that controls the maximum number of concurrent objects. The
  default value is five.
  <br>
  (Bruno Turcksin, 2020/05/29)
 </li>

 <li>
  New: Add multigrid transfer operators for distributed polynomial and
  global coarsening.
  <br>
  (Peter Munch, Laura Prieto Saavedra, 2020/05/29)
 </li>

 <li>
  Improved: step-28 now uses tasks instead of threads.
  <br>
  (David Wells, 2020/05/28)
 </li>

 <li>
  New: The class Particles::DataOut can now output particle properties as
  scalars, vectors, or tensors, depending on the arguments handed over to the
  Particles::DataOut::build_patches() function.
  <br>
  (Rene Gassmoeller, 2020/05/27)
 </li>

 <li>
  New: When executing a task on a separate thread, if that task ends
  with throwing an exception instead of returning a value, then this
  exception will be obtained when you wait for the task using
  Threads::Task::join() or Threads::Task::return_value().
  <br>
  (Wolfgang Bangerth, 2020/05/27)
 </li>

 <li>
  New: GridTools::Cache::get_locally_owned_cell_bounding_boxes_rtree() extracts
  a tree of bounding boxes covering the locally owned cells of a triangulation.
  This can be used in conjunction with GridTools::Cache::get_covering_rtree() to
  make approximate geometrical queries on who owns what spatial region.
  <br>
  (Luca Heltai, 2020/05/26)
 </li>

 <li>
  New: pack_rtree_of_indices() and IndexableGetterFromIndices allow to construct an
  RTree object that stores indices to existing containers of indexable objects.
  <br>
  (Luca Heltai, 2020/05/24)
 </li>

 <li>
  New: The class ParticleHandler now provides a signal 'signals.particle_lost'
  that is triggered whenever a particles can not be associated with a cell while
  calling its function sort_particles_into_subdomains_and_cells().
  <br>
  (Rene Gassmoeller, 2020/05/25)
 </li>

 <li>
  Bugfix: hp::Refinement::choose_p_over_h() now works in parallel.
  <br>
  (Marc Fehling, 2020/05/22)
 </li>

 <li>
  New: There is now a second overload for
  Particles::ParticleAccessor::set_properties() that takes an ArrayView
  as argument.
  <br>
  (Wolfgang Bangerth, 2020/05/22)
 </li>

 <li>
  Removed: All headers under <code>base/std_cxx11/</code> have been removed.
  <br>
  (David Wells, 2020/05/21)
 </li>

 <li>
  Changed: In many other classes that are templated on `dim` and
  `spacedim`, the second template argument `spacedim` had a default
  value equal to `dim`. Particles::DataOut did not, but now does.
  <br>
  (Wolfgang Bangerth, 2020/05/21)
 </li>

 <li>
  New: A new BoundingBoxDataOut class is available, to output a collection
  of objects that can be converted to BoundingBox objects in graphical format.
  <br>
  (Luca Heltai, 2020/05/19)
 </li>

</ol>

*/
