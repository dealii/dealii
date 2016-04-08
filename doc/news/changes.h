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
 <li> Fixed: Meshworker::Assembler::ResidualSimple now also works for
 multiple blocks if no constraints are given.
 <br>
 (Daniel Arndt, 2016/04/08)
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
