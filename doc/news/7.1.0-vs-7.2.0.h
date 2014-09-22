// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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
 * @page changes_between_7_1_and_7_2 Changes between Version 7.1 and 7.2

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
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
<li> Changed/fixed: Several operations on Trilinos vectors that change the
elements of these vectors were allowed accidentally for vectors that have
ghost elements. This is a source of errors because a change on one
MPI process will not show up on a different processor. Consequently, we
intended to disallow all functions that modify vectors with ghost elements
but this was not enforced for all of these functions. This is now fixed
but it may lead to errors if your code relied on the existing behavior. The
way to work around this is to only ever modify fully distributed vectors,
and then copy it into a vector with ghost elements.
<br>
(Wolfgang Bangerth, 2012/08/06)

<li> Changed: The argument material_id of the estimate-functions of
the KellyErrorEstimator class is now of type types::material_id
with default value numbers::invalid_material_id, instead of type
unsigned int with default value numbers::invalid_unsigned_int. This
should not make a difference to most users unless you specified the
argument's default value by hand.
<br>
(Wolfgang Bangerth, Christian Goll 2012/02/27)

<li>
The member functions Triangulation::set_boundary and
Triangulation::get_boundary now take a types::boundary_id instead of
an unsigned int as argument. This now matches the actual data type
used to store boundary indicators internally.
<br>
(Wolfgang Bangerth, Christian Goll 2012/02/27)
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
<li>
New: We now support parallel output using HDF5/xdmf.
<br>
(Eric Heien, Timo Heister, 2012/08/28)

<li>
New: We are now compatible with Trilinos 10.4.2, 10.8.5, and 10.12.2. See the
readme for more information.
<br>
(Timo Heister, 2012/08/24)

<li>
Changed: To make the naming of types defined in namespace types
consistent, we have renamed types::subdomain_id_t to
types::subdomain_id. The old name has been retained as a typedef
but is now deprecated.
<br>
(Wolfgang Bangerth, 2012/08/22)

<li>
Changed: Unify the concept of compress() for all linear algebra
objects. Introduce type VectorOperation to decide between
add and insert. Implement also for serial vectors. Note:
this breaks distributed::vector::compress(bool). See
@ref GlossCompress for more information.
<br>
(Timo Heister, 2012/08/13)

<li>
Changed: Support for the METIS 4.x has been replaced with support for
METIS 5.x. Use <code>--with-metis=path/to/metis</code> to configure
with METIS 5.x.
<br>
(Stefano Zampini, Toby D. Young 2012/08/13)

<li>
Changed: numerics/vectors.h is now called numerics/vector_tools.h and
numerics/matrices.h is now called numerics/matrix_tools.h The old files are
deprecated.
<br>
(Timo Heister 2012/08/09)

<li>
New: officially added support for clang 3.1 or newer.
<br>
(Timo Heister and Wolfgang Bangerth, 2012/08/07)

<li>
Changed: PETSc linking now prefers to use the libpetsc.so generated
by PETSc starting from version 3.1+. This fixes the problem
of linker errors on recent gcc/ubuntu versions.
<br>
(Timo Heister, 2012/08/07)

<li>
Fixed: On some 64-bit systems, we build deal.II with the <code>-m64</code>
flag but forgot to build UMFPACK with this flag as well, leading to
linker errors. This is now fixed.
<br>
(Wolfgang Bangerth, 2012/07/31)

Fixed: The Intel compiler, when using MPI, wants that <code>mpi.h</code>
is included before header files like <code>stdio.h</code>. This can't
be ensured in general because the inclusion might be indirectly, but
we now work around the problem in other ways.
<br>
(Timo Heister, Wolfgang Bangerth, Michael Thomadakis, 2012/07/26)

<li>
Fixed: On some systems, the p4est library we use for distributed
parallel computations installs its libraries into a <code>lib64/</code>
directory instead of the usual <code>lib/</code>. deal.II can now deal
with this.
<br>
(Wolfgang Bangerth, 2012/07/25)

<li>
New: step-43 is an extension of step-21 that shows efficient methods
to solve multi-phase flow.
<br>
(Chih-Che Chueh, Wolfgang Bangerth, 2012/06/06)

<li>
New: step-15 has been replaced by a program that demonstrates the
solution of nonlinear problem (the minimal surface equation) using
Newton's method.
<br>
(Sven Wetterauer, 2012/06/03)

<li>
New: step-48 demonstrates the solution of a nonlinear wave equation
with an explicit time stepping method. The usage of Gauss-Lobatto
elements gives diagonal mass matrices, which obviates the solution of
linear systems of equations. The nonlinear right hand side is
evaluated with the matrix-free framework.
<br>
(Katharina Kormann and Martin Kronbichler, 2012/05/05)

<li>
New: step-37 shows how the matrix-free framework can be utilized to
efficiently solve the Poisson equation without building global
matrices. It combines fast operator evaluation with a multigrid solver
based on polynomial Chebyshev smoother.
<br>
(Katharina Kormann and Martin Kronbichler, 2012/05/05)

<li>
New: A new matrix-free interface has been implemented. The framework
is parallelized with MPI, TBB, and explicit vectorization instructions
(new data type VectorizedArray). The class MatrixFree caches
cell-based data in an efficient way. Common operations can be
implemented using the FEEvaluation class.
<br>
(Katharina Kormann and Martin Kronbichler, 2012/05/05)

<li>
New: step-44 demonstrates one approach to modeling large
deformations of nearly-incompressible elastic solids. The
elastic response is governed by a non-linear, hyperelastic
free-energy function. The geometrical response is also
nonlinear, i.e., the program considers finite deformations.
<br>
(Andrew McBride and Jean-Paul Pelteret, 2012/04/25)

<li>
New: step-41 demonstrates solving the obstacle problem,
a variational inequality.
<br>
(Joerg Frohne, 2012/04/22)

<li>
Changed: The version of BOOST we ship with deal.II has been upgraded
to 1.49.0.
<br>
(Martin Kronbichler, 2012/04/07)

<li>
New: We have added a brief section to the step-12 tutorial programs that
compares the DG solution computed there with one that one would obtain by
using a continuous finite element.
<br>
(Wolfgang Bangerth, 2012/03/25)

<li>
New: Added support for codimension 2, i.e. for dim-dimensional objects
embedded into spacedim=dim+2 dimensional space.
<br>
(Sebastian Pauletti, 2012/03/02)

<li> Changed: Material and Boundary indicators have been both of the
type unsigned char. Throughout the library, we changed their datatype
to <code>types::material_id</code>
resp. <code>types::boundary_id</code>, both typedefs of unsigned
char. Internal faces are now characterized by
numbers::internal_face_boundary_id(<code>=static_cast@<types::boundary_id@>(-1)</code>)
instead of 255, so we get rid of that mysterious number in the source
code.  Material_ids are also assumed to lie in the range from 0 to
numbers::invalid_material_id-1 (where <code>numbers::invalid_material_id =
static_cast<types::material_id>(-1)</code>). With this change, it is now
much easier to extend the range of boundary or material ids, if
needed.
<br>
(Christian Goll 2012/02/27)

<li> New: Functions like FEValues::get_function_values have long been
able to extract values from pretty much any vector kind given to it (e.g.
of kind Vector, BlockVector, PETScWrappers::Vector, etc). The list of
allowed "vector" types now also includes IndexSet, which is interpreted
as a vector of elements that are either zero or one, depending on whether
an index is in the IndexSet or not.
<br>
As a consequence of this, the DataOut::add_data_vector functions now also
work for such types of vectors, a use demonstrated in step-41.
<br>
(Wolfgang Bangerth, 2012/02/14)

<li> New: It has long been a nuisance that the deal.II vector classes
could only be accessed using <code>operator()</code> whereas the C++
<code>std::vector</code> class required <code>operator[]</code>. This
diminished the usefulness of template code. Historically, the reason
was that the deal.II vector classes should use the same operator as
the matrix classes, and C++ does not allow to use <code>operator[]</code>
for matrices since this operator can only take a single argument.
<br>
In any case, all deal.II vector classes now support both kinds of
access operators interchangeably.
<br>
(Wolfgang Bangerth, 2012/02/12)

<li> Fixed: Linking shared libraries on PowerPC systems (e.g. on
BlueGene systems) failed due to a miscommunication between compiler
and linker. This is now worked around.
<br>
(Aron Ahmedia, Wolfgang Bangerth, 2012/02/06)

<li> New: There is now a distributed deal.II vector class
parallel::distributed::Vector that can be used with MPI. The
vector is based on a contiguous locally owned range and allows easy
access of ghost entries from other processors. The vector interface is
very similar to the non-distributed class Vector<Number>.
<br>
(Katharina Kormann, Martin Kronbichler, 2012/01/25)

<li> Fixed: The <code>common/scripts/make_dependencies</code> program
now behaves like the C++ compiler when
searching include paths for <code># include "..."</code> directives.
<br>
(Timo Heister, 2012/01/10)

<li> Fixed: The Intel compiler complains that it can't copy Trilinos vector
reference objects, preventing the compiling of step-32. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/09)

<li> Fixed: Intel ICC 12.1 gets into trouble with BOOST because BOOST
believes that the compiler supports C++0x but one then still has to
specify the corresponding flag on the command line to avoid compiler
errors. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/06)

<li> Fixed: On some systems, <code>mpiCC</code> turns out to alias the
C compiler, not the C++ compiler as expected. Consequently, try to use
<code>mpic++</code> or <code>mpicxx</code> before <code>mpiCC</code> as
these should really be unambiguous.
<br>
(Wolfgang Bangerth, 2011/11/05)

<li> Fixed: Intel's ICC compiler identifies itself as <code>icpc version
12.1.0 (gcc version 4.2.1 compatibility)</code> which we mistook as being
GCC version 4.2. The same is true for the Intel C compiler. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/05)

<li> Fixed: deal.II could not be compiled with gcc 4.6.1 when MPI is
enabled due to a missing include file in file
<code>source/base/utilities.cc</code>. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/25)
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
<li> Changed: To support Trilinos version 10.12.x we need to cache the
result of has_ghost_elements() in parallel vectors at construction time
of the vector. Starting from 10.12 Trilinos will communicate in this
call, which therefore only works if called from all CPUs.
<br>
(Timo Heister, 2012/08/22)

<li> New: The copy constructor of FullMatrix from IdentityMatrix
used to be explicit, but that didn't appear to be necessary in hindsight.
Consequently, it is now a regular copy constructor.
<br>
(Wolfgang Bangerth, 2012/08/14)

<li> New: The Patterns::Map pattern allows to describe maps from keys
to values in input files.
<br>
(Wolfgang Bangerth, 2012/08/01)

<li> Fixed: DoFTools::make_zero_boundary_constraints now also works for parallel distributed triangulations.
<br>
(Wolfgang Bangerth, 2012/07/24)

<li> Fixed: GridTools::find_active_cell_around_point() works now also if the cell in which the point we look for lies is not adjacent to the closest vertex of p. The tests bits/find_cell_8 and _9 illustrate this.
<br>
(Wolfgang Bangerth, Christian Goll 2012/07/20)

<li> Fixed: Using the SolutionTransfer class with hp::DoFHandler
and on meshes where some cells are associated with a FE_Nothing element
could result in an error. This is now fixed.
<br>
(Wolfgang Bangerth, 2012/06/29)

<li> Fixed: The MappingQ1::transform_real_to_unit_cell function as
well as the equivalent ones in derived classes sometimes get into
trouble if they are asked to compute the preimage of this point
in reference cell coordinates. This is because for points outside
the reference cell, the mapping from unit to real cell is not
necessarily invertible, and consequently the Newton iteration to
find the preimage did not always converge, leading to an exception.
While this is not entirely wrong (we could, after all, not compute
the desired quantity), not all callers of this function were prepared
to accept this result -- in particular, the function
CellAccessor<3>::point_inside should have really just returned false
in such cases but instead let the exception so generated propagate
through. This should now be fixed.
<br>
(Wolfgang Bangerth, Eric Heien, Sebastian Pauletti, 2012/06/27)

<li> Fixed: The function VectorTools::compute_no_normal_flux_constraints had
a bug that led to an exception whenever we were computing constraints for
vector fields located on edges shared between two faces of a 3d cell if those
faces were not parallel to the axes of the coordinate system. This is now fixed.
<br>
(Wolfgang Bangerth, Jennifer Worthen, 2012/06/27)

<li>
Fixed: Due to an apparent bug in autoconf, it was not possible to
override the <code>F77</code> environment variable to select anything
else than gfortran. This is now fixed.
<br>
(Wolfgang Bangerth, 2012/06/08)

<li> Fixed: TrilinosWrappers::VectorBase::swap() is now working as expected. (thanks Uwe K&ouml;cher)
<br>
(Timo Heister 2012/07/03)

<li> Fixed: Some instantiations for
DerivativeApproximation::approximate_derivative_tensor() were missing.
<br>
(Timo Heister 2012/06/07)

<li> New: The finite element type FE_DGQArbitraryNodes is now
working also in codimension one spaces.
<br>
(Luca Heltai, Andrea Mola 2012/06/06)

<li> Fixed: Computing the $W^{1,\infty}$ norm and seminorm in
VectorTools::integrate_difference was not implemented. This is now
fixed.
<br>
(Wolfgang Bangerth 2012/06/02)

<li> Fixed: A problem in MappingQ::transform_real_to_unit_cell
that sometimes led the algorithm in this function to abort.
<br>
(Wolfgang Bangerth 2012/05/30)

<li> New: The function DataOutInterface::write_pvd_record can be used
to provide Paraview with metadata that describes which time in a
simulation a particular output file corresponds to.
<br>
(Marco Engelhard 2012/05/30)

<li> Fixed: A bug in 3d with hanging nodes in GridTools::find_cells_adjacent_to_vertex()
that caused find_active_cell_around_point() to fail in those cases.
<br>
(Timo Heister, Wolfgang Bangerth 2012/05/30)

<li> New: DoFTools::make_periodicity_constraints implemented which
inserts algebraic constraints due to periodic boundary conditions
into a ConstraintMatrix.
<br>
(Matthias Maier, 2012/05/22)

<li> New: The GridIn::read_unv function can now read meshes generated
by the Salome framework, see http://www.salome-platform.org/ .
<br>
(Valentin Zingan, 2012/04/27)

<li> New: There is now a second DoFTools::map_dofs_to_support_points
function that also works for parallel::distributed::Triangulation
triangulations.
<br>
(Wolfgang Bangerth, 2012/04/26)

<li> New: There is now a second DoFTools::extract_boundary_dofs
function that also works for parallel::distributed::Triangulation
triangulations.
<br>
(Wolfgang Bangerth, 2012/04/26)

<li> New: The FullMatrix::extract_submatrix_from, FullMatrix::scatter_matrix_to,
FullMatrix::set functions are new.
<br>
(Andrew McBride, Jean Paul Pelteret, Wolfgang Bangerth, 2012/04/24)

<li> Fixed:
The method FEValues<dim>::inverse_jacobian() previously returned the transpose of the
inverse Jacobians instead of just the inverse Jacobian as documented. This is now fixed.
<br>
(Sebastian Pauletti, Katharina Kormann, Martin Kronbichler, 2012/03/11)

<li> Extended:
SolutionTransfer is now also able to transfer solutions between hp::DoFHandler where
the finite element index changes during refinement.
<br>
(Katharina Kormann, Martin Kronbichler, 2012/03/10)

<li> Changed:
A new method to determine an initial guess for the Newton method was coded
in MappingQ::transform_real_to_unit_cell.
The code in transform_real_to_unit_cell was cleaned a little bit and a new code
for the @<2,3@> case was added.
<br>
(Sebastian Pauletti, 2012/03/02)

<li> Changed:
In the context of codim@>0, Mapping::transform would require different inputs
depending on the mapping type.
For mapping_covariant, mapping_contravariant the input is DerivativeForm<1, dim, spacedim>
but for mapping_covariant_gradient,  mapping_contravariant_gradient the input is Tensor<2,dim>.
<br>
(Sebastian Pauletti,  2012/03/02)

<li> New:
A new class DerivativeForm was added.
This class is supposed to represent the derivatives of a mapping.
<br>
(Sebastian Pauletti, 2012/03/02)

<li> Fixed: TrilinosWrappers::Vector::all_zero() in parallel.
<br>
(Timo Heister, J&ouml;rg Frohne, 2012/03/06)

<li> New: GridGenerator::quarter_hyper_shell() in 3d.
<br>
(Thomas Geenen, 2012/03/05)

<li> New: DataOut::write_vtu_in_parallel(). This routine uses MPI I/O to write
a big vtu file in parallel.
<br>
(Timo Heister, 2012/02/29)

<li> Fixed: parallel::distributed::Triangulation::clear() forgot
to update the number cache of this class, leading to wrong results
when calling functions like
parallel::distributed::Triangulation::n_global_active_cells();
<br>
(Wolfgang Bangerth, 2012/02/20)

<li> Improved: FEFieldFunction allows now for the computation of Laplacians.
<br>
(Christian Goll, 2012/02/16)

<li> New: The function IndexSet::fill_binary_vector creates a numerical
representation of an IndexSet containing zeros and ones.
<br>
(Wolfgang Bangerth, 2012/02/12)

<li> New: The function IndexSet::clear resets an index set to be empty.
<br>
(Wolfgang Bangerth, 2012/02/12)

<li> New: There are now global functions l1_norm() and linfty_norm() that compute
the respective norms of a rank-2 tensor.
<br>
(Wolfgang Bangerth, 2012/02/08)

<li> New: VectorTools::interpolate_to_different_mesh implemented which interpolates between
     DoFHandlers with different triangulations based on a common coarse grid.
<br>
(Matthias Maier, 2012/02/08)

<li> Improved: DoFTools::map_dofs_to_support_points() now also works within the hp framework.
<br>
(Christian Goll, 2012/02/02)

<li> Improved: There is now a constructor for FESystem that allows to
create collections of finite elements of arbitrary length.
<br>
(Jason Sheldon, 2012/01/27)

<li> Improved: VectorTools::point_value() now also works within the hp framework.
<br>
(Christian Goll, 2012/01/26)

<li> Fixed: GridTools::find_active_cell_around_point() for the hp-case works now also with MappingCollections containing only one mapping, as is the standard case in other functions using hp.
<br>
(Christian Goll, 2012/01/26)

<li> Fixed: parallel::distributed::refine_and_coarsen_fixed_fraction()
contained a rounding bug that often produced wrong results.
<br>
(Timo Heister, 2012/01/24)

<li> Improved: Utilities::break_text_into_lines now also splits the string at '\\n'.
<br>
(Timo Heister, 2012/01/17)

<li> Fixed: When <code>./configure</code> does not detect the presence
of <code>zlib</code>, writing output in VTU format failed to produce
a valid output file.
<br>
(Timo Heister, 2012/01/03)

<li> Improved: <code>PETScWrappers::SolverXXX</code> class was
restricted to using default solver options for the KSP only. It is now
possible to override those by using PETSc command-line options
<code>-ksp_*</code>; giving greater flexibility in controlling PETSc
solvers. (See the class's documentation).
<br>
(Vijay S. Mahadevan, 2011/12/22)

<li> New: The GridIn class now also reads the GMSH format 2.2 as written by
GMSH 2.5.
<br>
(Vijay S. Mahadevan, Wolfgang Bangerth, 2011/12/19)

<li> Improved: The GridRefinement::refine_and_coarsen_optimize function
assumed that the expected convergence order was 2. It has now gotten an
argument by which the user can prescribe a different value. A bug has also
been fixed in which the function incorrectly assumed in its algorithm that
the mesh was two-dimensional.
<br>
(Christian Goll, 2011/12/16)

<li> Fixed: Restriction and prolongation didn't work for elements of
kind FE_Nothing. Consequently, many other parts of the library also
didn't work, such as the SolutionTransfer class. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/12/15)

<li> Fixed: The DerivativeApproximation class now works for distributed computations.
<br>
(Timo Heister, 2011/12/15)

<li> Changed: The ExcMessage exception class took an argument of
type <code>char*</code> that was displayed when the exception
was raised. However, character pointers are awkward to work with
because (i) they can not easily be composed to contain things like
file names that are only known at run-time, and (ii) the string
pointed to by the pointer had to live longer than the local expression
in which the exception is generated when using the AssertThrow macro
(because, when we create an exception, the exception object is passed
up the call stack until we find a catch-clause; at that time, however,
the scope in which the exception object was created has long been left).
This restriction made it impossible to construct the message using std::string
and then just do something like <code>(std::string("file: ")+filename).c_str()</code>.
<br>
To remedy this flaw, the type of the argument to ExcMessage has been
changed to std::string since objects of this type are readily copyable
and therefore live long enough.
<br>
(Wolfgang Bangerth, 2011/12/14)

<li> New: Setting up a class derived from DataPostprocessor required some
pretty mechanical steps in which one has to overload four member functions.
For common cases where a postprocessor only computes a single scalar or
a single vector, this is tedious and unnecessary. For these cases, the
new classes DataPostprocessorScalar and DataPostprocessorVector provide
short cuts that make life simpler.
<br>
(Wolfgang Bangerth, 2011/12/14)

<li> Changed: The DataPostprocessor class previously required users of this
class to overload DataPostprocessor::get_names(),
DataPostprocessor::get_data_component_interpretation()
and DataPostprocessor::n_output_variables(). The latter function is redundant
since its output must equal the length of the arrays returned by the
first two of these functions. It has therefore been removed.
<br>
(Wolfgang Bangerth, 2011/12/14)

<li> Improved: Objects of the type LogStream::Prefix can now be used
as a safe implementation of the push and pop mechanism for log
prefices.
<br>
(Guido Kanschat, 2011/12/09)

<li> New: IndexSet::add_indices(IndexSet).
<br>
(Timo Heister, 2011/12/09)

<li> Fixed: Finite element Hessians get computed for codimension one,
at least for FE_Poly derived classes.
<br>
(Guido Kanschat, 2011/12/07)

<li> Changed: Output of ParameterHandler::print_parameters with argument
ParameterHandler::LaTeX was not particularly readable. The output has
therefore been rewritten to be more structured and readable.
<br>
(Wolfgang Bangerth, 2011/11/28)

<li> Fixed: The TimerOutput class set the alignment of output to right-aligned
under some circumstances, but didn't reset this to the previous value at the
end of output. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/11/28)

<li> New: The copy constructor of the SparseMatrixEZ function now works
(rather than throwing an exception) if the copied matrix has size zero.
This is in accordance to the other matrix classes.
<br>
(Wolfgang Bangerth, 2011/11/15)

<li> New: The class ScalarFunctionFromFunctionObject provides a quick way to
convert an arbitrary function into a function object that can be passed
to part of the library that require the Function@<dim@> interface.
The VectorFunctionFromScalarFunctionObject class does a similar thing.
<br>
(Wolfgang Bangerth, 2011/11/15)

<li> New: The DoFTools::count_dofs_per_block function now also works
for objects of type hp::DoFHandler.
<br>
(Jason Sheldon, 2011/11/13)

<li> New: FETools::get_fe_from_name() can now return objects of type FE_Nothing.
<br>
(Scott Miller, Jonathan Pitt, 2011/11/10)

<li> New: Implementation of an alternative handling of
inhomogeneous constraints in ConstraintMatrix. This is controlled with
a new parameter use_inhomogeneities_for_rhs in
distribute_local_to_global() and determines whether the correct or
zero values (this was the case before and still is the default) are
kept in the linear system during the solution process.
<br>
(J&ouml;rg Frohne, 2011/11/01)

<li> Fixed: SparseMatrix::mmult and SpareMatrix::Tmmult had a number of
issues that are now fixed: (i) rebuilding the sparsity pattern was allowed
even if several of the matrices involved in these operations shared a
sparsity pattern; (ii) the functions had a vector argument that had a default
value but the default value could not be used because it wasn't used in a
template context deducible by the compiler.
<br>
(Wolfgang Bangerth, 2011/10/30)

<li> New:
parallel::distributed::Triangulation::mesh_reconstruction_after_repartitioning
setting which is necessary for save()/load() to be deterministic. Otherwise
the matrix assembly is done in a different order depending on the order of old
refinements.
<br>
(Timo Heister, 2011/10/26)

<li> New: TriaAccessor<>::minimum_vertex_distance().
<br>
(Timo Heister, 2011/10/25)

<li> New: TableHandler::print_text now supports not only printing column
keys above their own column, but also in a separate header, to make it simpler
for external plotting programs to skip this line.
<br>
(Wolfgang Bangerth, 2011/10/22)

<li> Fixed: Trying to write a TableHandler object that is empty resulted
in a segmentation fault. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/21)

<li> New: The TableHandler class can now pad columns that have only been
partially filled. See the documentation of the class for a description.
<br>
(Wolfgang Bangerth, 2011/10/18)

<li> Fixed: In TableHandler::print_text, it can happen that the function
wants to print an empty string as the element of the table to be printed.
This can confuse machine readers of this table, for example for visualization,
since they then do not see this column in that row. To prevent this, we now
print <code>""</code> in such places.
<br>
(Wolfgang Bangerth, 2011/10/18)

<li> Fixed: Using Trilinos versions 10.4 and later on Debian failed to
configure due to a different naming scheme of Trilinos libraries on
Debian. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/17)

<li> Changed: The TableHandler class has been changed significantly
internally.  It could previously store arbitrary values (though in practice,
only int, unsigned int, double and std::string were implemented). The class is
now restricted to this particular set of types. On the other hand, the
TableHandler class can now be serialized.
<br>
(Wolfgang Bangerth, 2011/10/17)

<li> Fixed: searching in the doxygen documentation.
<br>
(Timo Heister, 2011/10/13)

<li> New: parallel::distributed::Triangulation::save()/load() to store
the refinement information to disk. Also supports saving solution vectors
using the SolutionTransfer class.
<br>
(Timo Heister, 2011/10/12)

<li> Fixed: The DataOut_DoFData::merge_patches did not compile with newer compilers.
This is now fixed.
<br>
(Wolfgang Bangerth, 2011/10/11)
</ol>


*/
