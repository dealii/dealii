// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

/**
@page changes_after_8_4_1 Changes after Version 8.4.1

<p>
This is the list of changes made after the release of deal.II version
8.4.1. All entries are signed with the names of the authors.
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
  <li> Improved: MeshWorker treats periodic faces as interior faces.
  <br>
  (Daniel Arndt, 2016/04/18)
  </li>

  <li> Changed: FlatManifold takes as argument a periodicity option. This
  used to be a Point<dim>, but it should have been a Tensor<1,dim>. This
  is now changed.
  <br>
  (Luca Heltai, 2016/04/09)
  </li>

  <li> Changed: The default nodal point distribution of FE_Q, FE_DGQ,
  FE_Q_DG0, FE_Q_Bubbles, and FE_TraceQ has been changed from equidistant
  points to the node points of the corresponding Gauss-Lobatto quadrature
  formula. For degrees one and two, the Gauss-Lobatto quadrature is
  equidistant and thus the unit support points are as before. However, the
  Gauss-Lobatto points are more dense towards the element boundaries at higher
  degrees. This gives well-conditioned interpolation at arbitrary orders and
  much more stable computations. While these node distribution was available
  before, it was not very visible and often lead to misunderstandings by
  inexperienced users. Most codes will not be affected by this change, even
  those using cubic and higher degree polynomials, apart from slightly
  different (better) interpolation behavior and different entries in solution
  vectors. If you explicitly need equidistant points, use the constructors
  <tt>FE_Q<dim>(QIterated<1>(QTrapez<1>(),degree))</tt> or
  <tt>FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapez<1>(),degree))</tt>.
  <br>
  (Martin Kronbichler, 2016/04/05)
  </li>

  <li> Removed: Support for the legacy <code>Make.global_options</code>
  file has been removed.
  <br>
  (Matthias Maier, 2016/03/17)
  </li>

  <li> Removed: Functions with names containing <code>boundary_indicator</code>
  have been removed. They had previously already been deprecated, and replaced
  by functions containing the string <code>boundary_id</code> instead, to keep
  with the style used for <code>material_id</code>, <code>subdomain_id</code>,
  etc.
  <br>
  (Wolfgang Bangerth, 2016/02/28)
  </li>

  <li> Changed: Many functions in VectorTools and MatrixTools now require
  matching data types between vectors, matrices, and Function arguments.
  <br>
  (Denis Davydov, 2016/02/27)
  </li>

  <li> Changed: ConstraintMatrix::distribute_local_to_global() and numerous
  functions in VectorTools namespace now require matching data types.
  This is done to correctly handle complex-valued case.
  <br>
  (Denis Davydov, 2016/02/22)
  </li>
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>

 <li> New: Added move operations to BlockIndices, BlockVectorBase and
 BlockVector; Vector move operations nullify old object instead of
 using swap.
 <br>
 (Daniel Shapero, 2016/04/13)
 </li>

 <li> New: Added TensorProductManifold to create new manifolds from two
 ChartManifold objects. This can be used, for example, to combine a
 2d manifold with a flat manifold for an extruded mesh.
 <br>
 (Timo Heister, 2016/04/12)
 </li>

 <li> New: Added New option in the read_ucd function of the GridIn class.
      A flag can now be assigned to the function, to decide wether the
      indicators specified in a UCD file should be interpreted as
      boundary_ids or as manifold_ids. This is particularly useful
      when the indicators refer to internal faces, for which
      boundary_ids cannot be used.
 <br>
 (Andrea Mola, 2016/04/11)
 </li>

 <li> New: Manifold objects were previously only used to compute the
 locations of individual new points on a manifold. Now, they are also
 used to compute tangent vectors (via Manifold::get_tangent_vector()), and this
 functionality provides the basis for computing normal vectors to manifolds
 as well.
 <br>
 In many cases, tangent vectors can be computed quite easily if the
 manifold has a functional description, i.e., if it can be
 represented via the ChartManifold class. In those cases, it is only
 necessary to overload the ChartManifold::push_forward_gradient()
 function that computes the derivatives of the push forward operation.
 <br>
 (Luca Heltai, Wolfgang Bangerth, 2016/04/08)
 </li>

 <li> New: Added CompositionManifold to create new manifolds from two
 ChartManifold objects. This can be used, for example, to rotate a
 cylindrical Manifold, or to make a cylinders with parabolic sides.
 <br>
 (Luca Heltai, 2016/04/09)
 </li>

 <li> New: Added a new Mapping class, MappingManifold, to use exact
 geometrical information extracted from the Manifold description instead
 of a polynomial approximation when computing transformations from the
 reference to the real cell. This class allows the computation of
 quadrature points, tangent vectors, and normal vectors which are exact
 with respect to the geometrical description, and it uses the underlying
 Manifold objects of the Triangulation. MappingManifold coincides with
 MappingQ1 for the FlatManifold descriptor.
 <br>
 (Luca Heltai, 2016/04/09)
 </li>

 <li> New: Added GridGenerator::torus() to generate the volume mesh of a
 torus in three dimensions and a manifold description TorusManifold to
 go with it.
 <br>
 (Timo Heister, 2016/03/21)
 </li>

 <li> New: Added GridTools::rotate() in three space dimensions.
 <br>
 (Timo Heister, 2016/03/18)
 </li>

 <li> New: Added unit tests for complex-valued PETSc and SLEPc.
 <br>
 (Toby D. Young, Denis Davydov, 2016/03/11)
 </li>

 <li> New: Added another scaling factor to Kelly error estimator, namely h_K.
 <br>
 (Denis Davydov, 2016/03/05)
 </li>

 <li> New: Added indent target to indent all headers and source
 files. Now you can do make (or ninja) indent inside the build
 directory.
 <br>
 (Alberto Sartori, 2016/03/02)
 </li>
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
 <li> Fixed: TrilinosWrappers::SparseMatrix will now exit early if there are no
 entries to add to the matrix. This usually occurs when zero elision is on. This
 fixes a bug where the matrix raises an exception if there are no entries to add
 to a matrix and the provided row and column values are not locally stored.
 <br>
 (David Wells, 2016/04/24)
 </li>

 <li> Fixed: TrilinosWrappers::MPI::Vector and TrilinosWrappers::Vector could
 access invalid memory in the reinit() method if the MPI communicator was
 deleted before termination of the program. This usually happened when using
 vectors from GrowingVectorMemory where a pool keeps vector alive. This has
 been fixed.
 <br>
 (Martin Kronbichler, 2016/04/23)
 </li>

 <li> Fixed: The methods TrilinosWrappers::SparseMatrix::(T)mmult previously
 produced invalid matrix sizes if the final matrix was non-square. This has
 been fixed.
 <br>
 (Martin Kronbichler, Daniel Jodlbauer, 2016/04/21)
 </li>

 <li> New: Added an optional string parameter to the ParameterHandler::read_input ()
 and ParameterHandler::read_input_from_string() functions.
 When a line which equals this string is encountered, the parsing of parameters
 is terminated.
 <br>
 (Denis Davydov, 2016/04/20)
 </li>

 <li> New: Added move operations to IndexSet.
 <br>
 (Daniel Shapero, 2016/04/19)
 </li>

 <li> Improved: The parallel loops in the deal.II Vector class for
 vector-vector operations have been revised for performance. This includes
 adjusting the minimum parallel grain size to 4096 vector entries and using an
 affinity partitioner provided by Threading Building Blocks for better data
 locality, especially on multi-socket systems.
 <br>
 (Martin Kronbichler, 2016/04/14)
 </li>

 <li> New: added ReinitHelper for PETSc. This is required by LinearOperator
 class to reinit vectors.
 <br>
 (Mauro Bardelloni, 2016/04/13)
 </li>

 <li> Fixed and improved: Fix algorithm for incomplete assignment of level
 subdomain ids for parallel geometric multigrid. Also optimize algorithms
 used for assignment and DoF communication.
 <br>
 (Timo Heister, Martin Kronbichler, 2016/04/12)
 </li>

 <li> Improved: DoFRenumbering::compute_Cuthill_McKee when used with
 distributed triangulations contained parts that scaled as the global problem
 size, rather than the processor-local size. This prevented its use with more
 than a few hundred cores when hanging node constraints were activated. This
 has been fixed.
 <br>
 (Martin Kronbichler, 2016/04/11)
 </li>

 <li> New: added hessenberg_signal and krylov_space_signal to SolverGMRES.
 These signals allow to retrieve the Hessenberg matrix and the basis vectors
 generated by the Arnoldi algorithm.
 <br>
 (Giuseppe Pitton, Luca Heltai, 2016/04/11)
 </li>

 <li> Fixed: Meshworker::Assembler::ResidualSimple now also works for
 multiple blocks if no constraints are given.
 <br>
 (Daniel Arndt, 2016/04/08)
 </li>

 <li> New: A move constructor has been added to Quadrature.
 <br>
 (Daniel Shapero, 2016/04/08)
 </li>

 <li> Fixed: The multigrid transfer performed invalid data accesses on
 multigrid hierarchies that define the coarse level as a level larger than
 0. This has been fixed.
 <br>
 (Martin Kronbichler, 2016/04/03)
 </li>

 <li> New: Add GridTools::remove_hanging_nodes() and
 GridTools::remove_anisotropy() in GridTools. GridTools::remove_hanging_nodes()
 detects cells with hanging nodes and refines the neighbours in the direction
 that removes hanging nodes or in every directions.
 GridTools::remove_anisotropy() refines a mesh until the resulting mesh is
 composed by cells with ratio between the extension in each coordinate
 direction lower than a fixed value.
 <br>
 (Mauro Bardelloni, 2016/03/28)
 </li>

 <li> New: When using C++11, a move constructor and assignment operator has
 been added to SparseMatrix, so that these objects can be returned from
 functions and packed into pairs and tuples.
 <br>
 (Daniel Shapero, 2016/03/27)
 </li>

 <li> New: The product of a rank-1 tensor (a vector) and a rank-2
 symmetric tensor (a symmetric matrix) is now defined and yields
 a rank-1 tensor (a vector). The opposite product was previously
 already defined.
 <br>
 (Wolfgang Bangerth, 2016/03/25)
 </li>

 <li> New: Triangulation::add_periodicity allows for accessing neighbors across
 periodic boundaries via new functions in TriaAccessor.
 <br>
 (Daniel Arndt, Ali Samii, 2016/03/23)
 </li>

 <li> Fixed: DoFHandler::locally_owned_dofs() could create a segmentation
 fault in cases where some processors do not own any cells. This was caused
 by an incorrect computation in DoFTools::locally_owned_dofs_per_subdomain().
 <br>
 (Wolfgang Bangerth, 2016/03/20)
 </li>

 <li> Improved: The distribution of degrees of freedom on multigrid levels,
 DoFHandler::distribute_mg_dofs(), contained a few steps that scaled
 quadratically in the number of local cells for certain configurations. These
 steps have been replaced by linear complexity calls.
 <br>
 (Martin Kronbichler, 2016/03/18)
 </li>

 <li> New: Added custom target "relocate" to Mac OS X builds, that runs
 a script to make all paths absolute in the shared libraries included
 in the deal.II package (only enabled when building a package, and when
 including external libraries to the package)
 <br>
 (Luca Heltai, 2016/03/14)
 </li>

 <li> New: Add NURBSPatchManifold. This class is a child of ChartManifold and
 implements a manifold descriptor for the face of a CAD imported usign
 OpenCASCADE.
 <br>
 (Mauro Bardelloni, 2016/03/09)
 </li>

 <li> New: When using C++11, there is now a function Threads::new_task()
 that can take as an argument either a lambda function, or the result
 of a std::bind expression, or anything else that can be called as in a
 function call. There is also a similar function Threads::new_thread()
 that takes the same kind of argument.
 <br>
 (Wolfgang Bangerth, 2016/03/07)
 </li>

 <li> New: When using C++11, the function filter_iterators() allows to filter a
 range of iterators using predicates (like the ones defined in IteratorFilter).
 <br>
 (Bruno Turcksin, 2016/03/04)
 </li>

 <li> Fixed: The OpenCASCADE::push_forward_and_differential_forms()
 function is now able to change the direction of the normal vector
 according to Orientation() method.
 <br>
 (Mauro Bardelloni, 2016/03/02)
 </li>

 <li> Fixed: The function IndexSet::make_trilinos_map() now works if some
 processors have a contiguous range of indices and others do not.
 <br>
 (Bruno Turcksin, 2016/02/17)
 </li>

 <li> Updated: step-44 has been been expressed in a more dimension independent
 manner, and can be now run in both 2-d and 3-d.
 <br>
 (Jean-Paul Pelteret, 2016/02/17)
 </li>

 <li> Fixed: FE_Nedelec elements up to polynomial order 12 can now be
 constructed.
 <br>
 (Jean-Paul Pelteret, 2016/02/12)
 </li>

 <li> Fixed: The GridTools::build_triangulation_from_patches() function now
 also copies the locations of vertices from the cells of the source
 triangulation to the triangulation that is built from the list of patch cells.
 <br>
 (Spencer Patty, 2016/02/11)
 </li>
</ol>

*/
