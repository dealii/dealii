// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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
 * @page changes_between_6_3_and_7_0 Changes between Version 6.3 and 7.0

<p>
This is the list of changes made between the deal.II releases listed above.
made to the three sub-libraries <a href="#base">base</a>,
<a href="#lac">lac</a>, and <a href="#deal.II">deal.II</a>, as well as
changes to the <a href="#general">general infrastructure,
documentation, etc</a>.
</p>

<p>
All entries are signed with the names of the author. Regular
contributor's names are abbreviated by WB (Wolfgang Bangerth), GK
(Guido Kanschat), RH (Ralf Hartmann).
</p>


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

<li> The public member variable Quadrature::n_quadrature_points has
been removed after being deprecated for nearly three
years. Quadrature.size() replaces its function, and the new version
more naturally preserves consistency.
<br>
(GK 2010/11/11)
</li>

<li>
        deal.II has had a somewhat quirky directory and library structure for
        historical reasons, but this has now changed: All include and source
        files are now under the top-level <code>include/</code>
        and <code>source</code> directories. Furthermore, we no longer build a
        plethora of libraries but only <code>libdeal_II.g.so</code> (debug
        version) and <code>libdeal_II.so</code> (optimized version). In
	particular, we no longer build different versions of the library for
	different space dimensions.
	<br>
	As a consequence, if your makefile makes any assumption on the
	location of deal.II include files or the name of the deal.II library
	it will need to be changed. The sample Makefiles have been
	updated for this.
	<br>
	(WB 2010/10/25)
	</li>

<li> The <code>Polynomials::Lagrange::generate_complete_basis</code>
function has been renamed to Polynomials::generate_complete_Lagrange_basis .
The function was previously the only member of a class, and a static
one on top of that, which did not make much sense.
<br>
(WB 2010/10/22)</li>

<li> The <code>QGauss2, QGauss3, ..., QGauss7</code> classes &mdash; deprecated
for more than 6 years already &mdash; have finally been removed.
You should use code like <code>QGauss@<dim@>(4)</code> instead.
<br>
(WB 2010/09/19)</li>


<li> The FE_Nedelec class had previously implemented the lowest order
when the value 1 was passed to the constructor. This has now been
adjusted: the lowest order now results from passing 0, making this
consistent with the FE_RaviartThomas class and following the convention
used in by Brezzi and Raviart (though not in the original paper of
Nedelec).
<br>
(Markus B&uuml;rg 2010/09/14)</li>

<li> The fields DoFHandler::tria and DoFHandler::selected_fe are now
private instead of protected. Inheriting classes can only access them
through DoFHandler::get_tria() and DoFHandler::get_fe(), respectively.
<br>
(GK 2010/08/16)</li>

</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li><p>New: The new step-38 program shows how to discretize and solve
  partial differential equations posed on curved manifolds embedded in
  higher dimensional spaces.
  <br>
  (Andrea Bonito, M. Sebastian Pauletti, 2011/01/02)
  </p></li>

  <li><p>New: deal.II now has a Qt based graphical user interface to edit
  input parameter files for the ParameterHandler class. The executable is
  located in <code>lib/bin/dealii_parameter_gui</code>.
  <br>
  (Martin Steigemann, 2010/12/20)
  </p></li>

  <li><p>New: A significant number of classes, primarily in the
  <code>base/</code> subdirectory, now provide the necessary member functions
  for serialization through BOOST's serialization library.
  <br>
  (Pradeep Rao, WB, 2010/12/06)
  </p></li>

  <li><p>Fixed: A good number of problems associated with meshes that are
  embedded in a higher space dimension (i.e. where <code>dim @< spacedim</code>)
  are now fixed.
  <br>
  (Sebastian Pauletti, Andrea Bonito, Luca Heltai, WB, 2010/12/06)
  </p></li>

  <li><p>Updated: The version of the <a href="http://www.threadingbuildingblocks.org">Threading
  Building Blocks (TBB)</a> shipped with deal.II has been updated
  to release 3.0 Update 3 (Commercially aligned version).
  <br>
  (WB, 2010/11/18)
  </p></li>

  <li><p>New: While 2d and 3d could mostly be handled with the
  same code, dimension-independent programming was always a bit
  different because the faces of 1d cells (i.e. vertices) did not
  have appropriate iterators defined. This has now been changed:
  there is now a class TriaAccessor<0,1,spacedim> (and corresponding
  DoFAccessor<0,DH>) that allows to write things like
  @code
    cell->face(0)->boundary_indicator();
    cell->face(1)->at_boundary();
    cell->face(1)->get_dof_indices(...);
  @endcode
  even if the cell corresponds to a one-dimensional triangulation where
  this was not previously possible.
  <br>
  (WB, 2010/11/11)
  </p></li>

  <li>
  <p>New: After more than 2 years of work, we have merged a branch on which
  we have implemented the functionality necessary to distribute meshes on
  hundreds, thousands, or more processors. An overview of the framework
  in which this is implemented is presented in the @ref distributed module,
  the @ref distributed_paper, as well as in step-40.
  <br>
  (Timo Heister, Martin Kronbichler, Wolfgang Bangerth 2010/10/23)
  </p>

  <li>
  <p>New: Documentation of how to handle constraints on degrees of freedom
  has been centralized in a new documentation module on @ref constraints.
  <br>
  (WB 2010/09/16)
  </p>

  <li>
  <p>Fixed: When using Trilinos and deal.II both with static libraries,
  a linker error would occur whenever a program linked both the 2d and
  3d libraries of deal.II. This is now fixed.
  <br>
  (WB 2010/07/23)
  </p>

  <li>
  <p>Fixed: On all non-linux platforms, if static libraries were selected
  and <code>./configure</code> was instructed to use Trilinos (which also
  was compiled with static libraries), a failure would occur. This should
  now be fixed.
  <br>
  (WB 2010/07/23)
  </p>

  <li>
  <p>Fixed: The step-33 tutorial program can not be built with GCC versions
  4.5.x. There are in fact two problems, one that pertains to uses of
  <code>std::make_pair</code> that don't work any more with the upcoming
  C++ 1x standard that GCC 4.5.x already follows, and some in which the
  Trilinos package Sacado is incompatible with GCC 4.5.x, at least up to
  Trilinos version 10.4.0. While the latter problem can only be fixed in
  future Trilinos versions, at least the former problem is solved in step-33.
  <br>
  (WB 2010/07/18)
  </p>

  <li>
  <p>Fixed: GCC version 3.4.0 failed to compile the file
  <code>deal.II/source/numerics/matrices.cc</code> with
  an internal compiler error. This has
  now been worked around.
  <br>
  (WB 2010/07/15)
  </p>

  <li>
  <p>Fixed: A problem in the Makefiles caused error messages when
  building under CygWin.
  <br>
  (GK 2010/07/12)
  </p>

  <li>
  <p>Fixed: GCC version 3.3.x failed to compile the files
  <code>lac/include/lac/precondition_block.h</code>,
  <code>deal.II/source/multigrid/mg_dof_handler.cc</code> and
  <code>examples/step-34/step-34.cc</code>. These problems have
  now been worked around.
  <br>
  (WB 2010/07/12)
  </p>

  <li>
  <p>Fixed: Some older 3.x versions of GCC crashed compiling the functions in
  namespace DoFRenumbering::boost. There is now a configuration time test
  that checks that the compiler accepts the constructs in question. If the
  compiler does not, then these functions are disabled.
  <br>
  (WB 2010/07/01)
  </p>

  <li>
  <p>Fixed: Linking with more than one of the deal.II 1d, 2d, or 3d libraries
  when using static libraries did not work. This is now fixed. However, due to
  GCC bug <a href="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=10591"
  target="_top">10591</a>, GCC versions prior to and including 4.1.x will
  still not work. Working with shared libraries was not and is not affected
  by this problem.
  <br>
  (WB 2010/07/01)
  </p>

  <li>
  <p>Updated: The version of <a href="http://www.boost.org/">boost</a>
   included in the <code>contrib/</code> directory has been updated
   to 1.43.0.
  <br>
  (WB 2010/06/30)
  </p>

  <li>
  <p>Fixed: GCC version 4.0.1 had a bug that prevented it from compiling
  release 6.3.0 because it apparently had an infinite loop allocating
  memory when compiling <code>fe_values.cc</code> in optimized mode. This
  problem had been fixed in GCC 4.0.2, but some versions of Mac OS X still use
  this GCC version in their Xcode environment. In any case, the code in
  deal.II has been changed to avoid this problem.
  <br>
  (WB 2010/06/30)
  </p>

  <li>
  <p>Fixed: Configuring with an external BOOST version did not work when
  using shared libraries since the test ran in the wrong order with respect
  to another configure test. This is now fixed.
  <br>
  (Bradley Froehle 2010/06/29)
  </p>

  <li>
  <p>Updated: The conversion tool in <code>contrib/mesh_conversion</code> that
  can read CUBIT output and convert it into something that is readable by
  deal.II has been updated.
  <br>
  (Jean-Paul Pelteret 2010/06/28)
  </p>

  <li>
  <p>
  Fixed: deal.II release 6.3.0 did not compile with Trilinos versions 9.x and
  10.0. This is now fixed.
  <br>
  (Martin Kronbichler, WB 2010/06/28)
  </p>
</ol>



<a name="base"></a>
<h3>base</h3>


<ol>

<li> A static member function Point::unit_vector() has been introduced
to simplify the task of creating a vector in one coordinate direction
independent of the dimension.
<br>
(GK 2010/11/17)
</li>

<li><p> There are now Patterns::FileName and Patterns::DirectoryName classes
  that can be used to indicate that a given parameter is supposed to be
  a file or directory name.
  <br>
  (Martin Steigemann 2010/10/25)
  </p>

  <li><p> New: The ParameterHandler class is now built on the
  <a href="http://www.boost.org" target="_top">boost</a>
  <code>property_tree</code> library which provides a much better
  foundation for extensions. In particular, the description of parameters
  can now be exported in XML and JSON formats for processing with external
  programs, for example graphical user interfaces. As a consequence of
  the re-write, a bug in ParameterHandler::print_parameters_section
  was fixed when using the LaTeX output format.
  <br>
  (WB 2010/09/09)
  </p>

  <li><p> Fixed: The ParameterHandler::set() functions allowed to set values that
  did not satisfy the pattern given during declaration of the parameter. This
  is now fixed: the functions now throw an exception.
  <br>
  (WB 2010/09/09)
  </p>

  <li><p> Fixed: The Patterns::Integer and Patterns::Double classes did
  not properly check that a value given matched the pattern because they
  ignored text after an initial match. This led to <code>"3.141"</code> and
  <code>"3.141..,-RXYZ"</code> to be recognized as valid integers and
  double values, respectively. This is now fixed.
  <br>
  (WB 2010/09/08)
  </p>

  <li><p> Fixed: The computation of quadrature points in the QGaussLobatto
  class uses a Newton method that was wrongly implemented. While the
  results were correct (at least for moderate orders), it required more
  iterations than necessary. This is now fixed.
  <br>
  (Andrea Bonito 2010/08/12)
  </p>

  <li><p> Changed: The DataOutBase::write_vtu function now writes data
  as 64-bit values, rather than 32-bit values.
  <br>
  (Scott Miller 2010/08/5)
  </p>

  <li><p> New: MappingQ and MappingQEulerian now support order > 1 also in
  codimension one. Step-34 has been modified to show how this works.
  <br>
  (Luca Heltai 2010/07/23-27)
  </p>

  <li><p> New: QGaussOneOverR now has a new constructor for arbitrary quadrature
  points and not only the vertices of the reference cell.
  <br>
  (Luca Heltai 2010/07/21)
  </p>
</ol>


<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li><p>New: Added a new constructor for PETScWrappers::Vector that takes
  an existing PETSc Vec that was created by the user and is only wrapped for
  usage.
  <br>
  (TH, 2010/12/03)
  </p></li>

  <li><p>New: When Hypre is available as a sub-package of PETSc, the new class
  PETScWrappers::PreconditionBoomerAMG provides the BoomerAMG implementation
  of the algebraic multigrid method of the Hypre package as a preconditioner.
  This functionality is used, for example, in step-40.
  <br>
  (Timo Heister 2010/11/13)
  </p></li>

  <li><p>New: The class LAPACKFullMatrix now has functions to compute the
  singular value decomposition of a matrix and its inverse.
  <br>
  (GK 2010/11/7)
  </p></li>

  <li><p>New: The classes RelaxationBlockSOR and RelaxationBlockSSOR
  implement overlapping Schwarz relaxation methods. Additionally,
  their base class RelaxationBlock and the helper class BlockList have
  been added to the library.
  <br>
  (GK 2010/10/19)
  </p></li>

  <li><p>Improved: The ConstraintMatrix::merge function now takes a second
  argument that indicates what should happen if the two objects to be
  merged have constraints on the very same degree of freedom.
  <br>
  (WB 2010/09/18)
  </p></li>

  <li><p>Updated: Changes to the <code>SLEPcWrappers</code> to preserve
  compatibility with SLEPc version 3.1. Main new updated features are
  new solver classes: (i) Power; and (ii) Generalized Davidson.
  <br>
  (Toby D. Young 2010/08/04)
  </p></li>

  <li><p>Fixed: SparseMatrix::precondition_SSOR and
  SparseMatrixEZ::precondition_SSOR had a bug, which made the SSOR
  method inaccurate for relaxation parameters different from unity.
  <br>
  (GK 2010/07/20)
  </p></li>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li><p>New: The GridTools::transform function used to work only on
  meshes that are not refined at all. It now also works on meshes
  that have been refined uniformly, though not on adaptively refined
  ones.
  <br>
  (WB, 2010/12/20)
  </p></li>
  
  <li><p>New: When computing errors using
  VectorTools::integrate_difference in codimension one problems, if
  you specified a norm that requires the computation of the gradients,
  you would get in trouble, because on codimension one manifolds we only
  have information on the <em>tangential</em> gradient. This is now
  fixed, by discarding the normal component of the provided function,
  before computing the difference with the finite element function.
  <br>
  (Luca Heltai, 2010/12/05)
  </p></li>
  
  <li><p>New: The VectorTools::interpolate_boundary_values function can now
  also be used to interpolate boundary values on meshes embedded in a higher
  dimensional space.
  <br>
  (WB, 2010/11/20)
  </p></li>

  <li><p>Fixed: Triangulation::execute_coarsening_and_refinement misplaced
  the new central vertex of a quad in a two-dimensional mesh embedded into
  higher dimensional space if said quad was at the boundary of the surface
  described by this triangulation. This is now fixed.
  <br>
  (WB, 2010/11/18)
  </p></li>

  <li><p>New: The GridTools::extract_boundary_mesh function can be used to
  generate a surface mesh from the boundary of a volume mesh.
  <br>
  (Sebastian Pauletti, 2010/11/15)
  </p></li>

  <li><p>Fixed: For DoF iterators, it was previously possible to write
  code like <code>*it1 = *it2</code>, presumably with the intent to
  copy the entire cell pointed to on the right hand side onto the cell
  pointed to at the left. However, this is not what happens since
  iterators are not pointers but only point to accessor classes. The
  assignment operator has therefore been removed.
  <br>
  (WB, 2010/11/12)
  </p></li>

  <li><p>New: The GridOut::write_gnuplot function now also works for meshes
  that are embedded in a higher
  dimensional space (i.e. if the codimension was greater than 1).
  <br>
  (WB, 2010/11/03)
  </p></li>

  <li><p>New: The class hp::FEFaceValues and hp::FESubfaceValues were not
  previously available for meshes that were embedded in a higher
  dimensional space (i.e. if the codimension was greater than 1). This is
  now fixed. As a consequence, the VectorTools::interpolate_boundary_values
  function is now also available for such meshes.
  <br>
  (WB, 2010/11/02)
  </p></li>

  <li><p>Fixed: The FEValuesExtractors::Vector class did not work when the dimension
  of the domain was not equal to the dimension of the space in which it is
  embedded. This is now fixed.
  <br>
  (Sebastian Pauletti, WB, 2010/11/01)
  </p></li>

  <li><p>Fixed: The methods VectorTools::project_boundary_values_curl_conforming
  set the boundary values now also H(curl)-conforming, when hanging node
  constraints are present.
  <br>
  (Markus Buerg, 2010/10/11)
  </p></li>

  <li><p>Changed: The DataPostprocessor functions now take an additional
  argument that indicates the location of an evaluation point. For backward
  compatibility, the old function signature still exists so that applications
  that overload one of the existing functions continue to work.
  The old signature has been deprecated, however, and will be removed in a
  future version.
  <br>
  (Scott Miller 2010/10/08)
  </p></li>

  <li><p>Changed: FETools is now a namespace rather than a class with only
  static member functions.
  <br>
  (WB 2010/09/25)
  </p></li>

  <li><p>New: FE_DGVector implements discontinuous elements based on
  vector valued polynomial spaces.
  <br>
  (GK 2010/09/17)
  </p></li>

  <li>
  <p>
  Fixed: The methods VectorTools::interpolate_boundary_values and
  VectorTools::project_boundary_values with ConstraintMatrix argument added
  inhomogeneities also to DoFs that were already constrained when entering the
  call without changing the list of DoFs that those were constrained to. This
  leads to wrong results for inhomogeneous constraints in case e.g.
  DoFTools::make_hanging_node_constraints was called before. Now, the correct
  way for combining hanging node constraints and Dirichlet conditions in one
  ConstraintMatrix is to first make the hanging node constraints, and then
  interpolate/project the boundary values. The latter operation will then
  simply not add any new constraints if a DoF is already constrained. The
  same is also true for VectorTools::compute_no_normal_flux_constraints.
  <br>
  (Luca Heltai, Martin Kronbichler, WB 2010/09/16)
  </p>

  <li>
  <p>
  Fixed: The method FEValuesViews::Vector::curl aborted the program in 2d under certain
  circumstances. This is now fixed.
  <br>
  (Andrea Bonito, WB 2010/09/14)
  </p>

  <li><p>Improved: MeshWorker::loop() now uses WorkStream for multithreading
  <br>
  (GK 2010/09/10)
  </p>

  <li>
  <p>
  Fixed: The method FEValuesViews::SymmetricTensor< 2, dim, spacedim >::get_function_divergences
  incorrectly calculated the divergence.
  <br>
  (Andrew McBride 2010/09/09)
  </p>

  <li>
  <p>
  New: Reimplementation of the parent function in the various DoFCellAccessors so
  that one can call the parent function from a DoFIterator and also gets a
  DoFIterator returned.
  <br>
  (Markus Buerg 2010/08/30)
  </p>

  <li>
  <p>
  New: FETools::compute_embedding_matrices now computes the embedding matrix
  for all refinement cases in parallel.
  <br>
  (Markus Buerg 2010/08/26)
  </p>

  <li>
  <p>
  New: There is now a function CellAccessor::parent so that one can do
  <code>cell-@>parent()</code>
  <br>
  (Markus Buerg 2010/08/26)
  </p>

  <li>
  <p>
  Improved: DoFHandler has a default constructor, so that it can be used in containers.
  <br>
  (GK 2010/08/16)
  </p>

  <li>
  <p>
  New: The functions VectorTools::project_boundary_values_curl_conforming
  are added. They can compute Dirichlet boundary conditions for Nedelec
  elements.
  <br>
  (Markus Buerg 2010/08/13)
  </p>

  <li>
  <p>
  Fixed: The function FEValuesViews::Vector::get_function_curls produced
  wrong results in some cases, because it erased the given vector first.
  This is now fixed.
  <br>
  (Markus Buerg 2010/08/13)
  </p>

  <li>
  <p>
  New: Ability to project second-order SymmetricTensor and first-order Tensor objects from the quadrature points to the support points of the cell using  FETools::compute_projection_from_quadrature_points
  <br>
  (Andrew McBride 2010/07/29)
  </p>

  <li>
  <p>
  Fixed: In some rather pathological cases, the function
  Triangulation::prepare_coarsening_and_refinement (which is called from
  Triangulation::execute_coarsening_and_refinement) could take very long
  if the flag Triangulation::eliminate_unrefined_islands was given in
  the mesh smoothing flags upon construction of the triangulation. This is
  now fixed.
  <br>
  (WB 2010/07/27)
  </p>

  <li><p> New: Brezzi-Douglas-Marini elements of arbitrary order in FE_BDM.
  <br>
  (GK 2010/07/19)
  </p>

  <li>
  <p>
  Fixed: The FEValues::get_cell() function was unusable from user code
  since its implementation used a class that was only forward declared
  and not visible at the point of instantiations in user code. This is now
  fixed.
  <br>
  (WB 2010/07/16)
  </p>

  <li>
  <p>New: The multilevel matrices for continuous elements
  can be built with the MeshWorker now.
  <br>
  (B&auml;rbel Janssen 2010/07/01)
  </p></li>

  <li>
  <p>
  Fixed: On some systems and compilers, the library could not be compiled
  because of a duplicate symbol in <code>MeshWorker::LocalResults</code>.
  This is now fixed.
  <br>
  (WB 2010/06/28)
  </p>

  <li>
  <p>
  Fixed: The output of the function
  FE_Q::adjust_quad_dof_index_for_face_orientation
  was wrong in 3d for polynomial orders of three or greater. This is now
  fixed.
  <br>
  (WB 2010/06/28)
  </p>
</ol>


*/
