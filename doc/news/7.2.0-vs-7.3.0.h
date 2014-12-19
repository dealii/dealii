// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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
 * @page changes_between_7_2_and_7_3 Changes between Version 7.2 and 7.3

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

<li> Changed: There are various changes to PETScWrappers::MPI::Vector
to minimize usage errors and to make the behavior similar to Trilinos:
objects with ghost elements are now read-only, the function
update_ghost_values() is no longer required (but called automatically
if needed). This requires some changes in user code.
<br>
(Timo Heister, 2013/02/16)

<li>Changed: Over the past few years, deal.II has accumulated a
number of things that we would like to change but that would introduce
incompatibility. Examples are inconsistent naming of functions or types,
or things put into the wrong namespace. In reality, while there are
about 100 examples of things we'd like to get rid of, most of them are
rather obscure to begin with and will not affect most user code.
Nevertheless, as a developer community, we are
very careful in making such incompatible changes.

All this said, it is sometimes necessary. We plan to create an incompatible
release 8.0 at a later time. To give users a chance to already see which
functions need to be changed, we have introduced markers into the deal.II
header files that identify which functions, classes or symbols are now
deprecated. If your compiler supports this, it will then warn if you are
using any of these. The documentation for each of these symbols discusses
the recommended alternative, which are in all cases already in place. In
other words, you can already change your code in such a way that it compiles
both with the current 7.3 release as well as with the future 8.0 release
in which these symbols will have been removed.
<br>
(Matthias Maier, Timo Heister, Wolfgang Bangerth, 2013/1/5)

<li>Changed: deal.II previously had two separate classes, DoFHandler
and MGDoFHandler, for standard and multilevel discretizations, respectively.
There were also corresponding accessor hierarchies for the cell and
face iterators. This turned out to be an obstacle for future developments.
Consequently, the functionality of the MGDoFHandler class has been merged
into the standard DoFHandler (which already was a base class of MGDoFHandler;
moving the functions down in the class hierarchy therefore does not
significantly affect how these classes are used). The MGDoFHandler class
still exists for now, but is deprecated and will be removed in deal.II 8.0.
Use the DoFHandler class instead which has obtained all the functionality of
the old MGDoFHandler class.

From a user perspective, there is only one significant change: The old
MGDoFHandler::distribute_dofs() function called DoFHandler::distribute_dofs()
and then distributed degrees of freedom on all levels of the multilevel
hierarchy. In the new scheme, this is now achieved using the function
DoFHandler::distribute_mg_dofs(). In other words, where you previously
used the class MGDoFHandler and called MGDoFHandler::distribute_dofs(), you
should now use the class DoFHandler and call first DoFHandler::distribute_dofs()
and then DoFHandler::distribute_mg_dofs().
<br>
(Markus B&uuml;rg, Timo Heister, Guido Kanschat, 2013/01/03)

<li>Removed: The Triangulation and DoFHandler classes had a great number
of functions that allowed to get the first and last iterators to cells,
faces, lines, quads and hexes overall and on each level of the mesh
hierarchy individually. While conceptually a nice idea to offer such a
rich set of functions, there are two factors that led us to drastically
reduce this set of functions: (i) A large interface makes it more
difficult to fully document what every function is doing and more
difficult to find what one is looking for. (ii) A large interface
makes it incredibly difficult to evolve the data structures on which
these functions operate. We felt that the many functions we have impeded
our ability to make changes to the internal storage of data. A survey
of the deal.II code base as well as some larger applications also shows
that most of these functions are rarely used, if at all.
<br>
Consequently, most of these functions have been removed. The only functions
that remain are the following (for each of the Triangulation and the
various DoFHandler classes):
<pre>
<code>
    cell_iterator        begin       (const unsigned int level = 0) const;
    active_cell_iterator begin_active(const unsigned int level = 0) const;
    cell_iterator        end         () const;
    cell_iterator        end         (const unsigned int level) const;
    active_cell_iterator end_active  (const unsigned int level) const;
</code>
</pre>
<br>
Codes that have previously used functions like <code>begin_active_line</code>
etc. to loop over the lines or quads of a triangulation need to be changed
to loop over cells, and on each cell loop over the lines or quads of this
cell. In most cases we have encountered (in deal.II or its testsuite) this
was a rather trivial modification. A case to watch out for is that the
old loop over all lines encountered all lines only once whereas one may
encounter it multiple times when looping over all cells and then the lines
of each cell. This case is easily avoided by flagging each treated line
using the @ref GlossUserFlags "user flags" associated with lines, quads,
and cells.
<br>
(Wolfgang Bangerth, 2012/09/22)

<li>New: In the past, deal.II used a std::vector of bools in many places
to denote component masks (see @ref GlossComponentMask) as well as for
block masks (see @ref GlossBlockMask). This was neither
descriptive (the data type does not indicate what it is supposed to
represent, nor whether the proper size for such an argument would be equal
to the number of degrees of freedom per cell, the number of vector components
of the finite element, or the number of degrees of freedom in total).
<br>
There are now new class types ComponentMask and BlockMask that are used in these places.
They are used both descriptively (as a return type of the function
FiniteElement::get_nonzero_components indicating the vector components within
which a given shape function is nonzero) as well as prescriptively (as
input arguments to functions such as those listed in the glossary entry
linked to above).
<br>
While the descriptive places are not backward compatible (they return a
ComponentMask which is not convertible to the std::vector of bools returned
before), most of the prescriptive places are in fact backward compatible
(because the std::vector of bool that was passed previously can
implicitly be converted to an object of type ComponentMask. The sole
exception is the function DoFTools::extract_dofs (and its multigrid
equivalent DoFTools::extract_level_dofs) that previously
could interpret its argument as either a component or
a block mask, depending on a boolean flag. This function now exists
in two different versions, one that takes a ComponentMask and
one that takes a BlockMask. Call sites need to be adjusted.
<br>
(Wolfgang Bangerth, 2012/09/22)

<li> Removed: the optional argument offset got removed from
DoFHandler and MGDoFHandler::distribute_dofs() because it was
never working correctly and it is not used.
<br>
(Timo Heister, 2012/09/03)
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>

<li> Changed: parallel sparse matrices now also require the use of
compress(VectorOperation) like vectors do. See the glossary for details.
Old functions are now deprecated.
<br>
(Timo Heister, 2013/02/25)

<li> New: step-49 demonstrates advanced techniques for mesh creation and
manipulation.
<br>
(Timo Heister, 2013/02/19)

<li> Fixed: Many places in the documentation have been made to consistently
use doxygen markup that indicates that this is code, so that doxygen will
cross-link the piece of code to class and function names. Many typos have also
been fixed.
<br>
(Felix Gruber, 2013/02/15)

<li> Fixed: Starting in release 7.1, we first built the deal.II shared libraries
in the local <code>/tmp</code> or similar directory rather than the final location
because linking becomes very slow over remotely mounted file systems. Unfortunately,
this schemes turns out not to work under Cygwin on Windows systems: executables
cannot link with libraries that have been moved/renamed after linking. If we are
running on Cygwin, we therefore revert to the old scheme.
<br>
(Wolfgang Bangerth, 2013/02/11)

<li> New: finite element FE_Q_DG0 that implements polynomials
of order k with an additional discontinuous constant function.
<br>
(Daniel Arndt, Timo Heister, 2013/01/07)

<li> step-6 now uses ConstraintMatrix::distribute_local_to_global()
instead of condense(), which is the preferred way to use a ConstraintMatrix
 (and the only sensible way in parallel).
<br>
(Timo Heister, 2012/11/02)

<li> Simplifications of the internal structures of Triangulation and
DoFHandler, in particular removal of specializations.
<br>
(Guido Kanschat, 2012/09/13)
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
<li> Fixed: The ArpackSolver interface to the ARPACK eigenvalue solver could
not be compiled with newer C++ compilers. This is now fixed.
<br>
(Juan Carlos Araujo Cabarcas, Wolfgang Bangerth, 2013/02/20)

<li> New: PETScWrappers::MPI::BlockVector now has a constructor and reinit
that takes a std::vector<IndexSet> (same interface as in Trilinos).
<br>
(Timo Heister, 2013/02/19)

<li> New: PETScWrappers::*Matrix::%add(other, factor) to
add a scaled other matrix to the current matrix.
<br>
(Jose Javier Munoz Criollo, 2013/02/19)

<li> New: GridGenerator::extrude_triangulation() allows
you to extrude a 2d mesh to turn it into a 3d mesh.
<br>
(Timo Heister, 2013/02/16)

<li> PETScWrappers::Vector::operator= and PETScWrappers::MPI::Vector::operator=
now call update_ghost_values()
automatically if necessary. This means that update_ghost_values()
does not need to be called from user code at all any more.
<br>
(Timo Heister, 2013/02/14)

<li> Fixed: VectorTools::interpolate did not work properly in 1d if
boundary indicators had been set to anything but the default (i.e.,
zero at the left and one at the right end of the domain). This was
a hold-over from the past when these were the only possible values.
This is now fixed.
<br>
(Kevin Dugan, Wolfgang Bangerth, 2013/02/14)

<li> Improved: The iterator class of the deal.II SparseMatrix class
and SparsityPattern have been revised for performance. Iterating over
a row of the matrix and querying the column index and the value is now
similarly fast as iterating over a vector.
<br>
(Martin Kronbichler, 2013/02/12)

<li> New: A new overload of BlockMatrixBase::add allows to add one block
matrix to another, with a scaling factor.
<br>
(Jean-Paul Pelteret, 2013/02/05)

<li> Fixed: The FEValues machinery silently accepted the case when the
mapped cell (or the cell geometry) were distorted. An assertion has been
added to the computation of the Jacobian determinants for the volume
element that aborts whenever the Jacobian determinant in a quadrature
point becomes too small or negative.
<br>
(Martin Kronbichler, 2013/01/18)

<li> Improved: SLEPcWrappers:: The interface to SLEPc has an improved
handle on SolverControl and solver data can now be extracted at run
time. An example usage has been added to step-36.
<br>
(Toby D. Young, 2013/01/18)

<li> Fixed: Various variants of the TrilinosWrappers::SparseMatrix::reinit
functions take a parameter <code>drop_tolerance</code> that allows to remove
small values from the matrix and replace them by zero instead. This was not
enforced for values on the diagonal of the matrix but only for off-diagonal
ones. This is now fixed.
<br>
(Wolfgang Bangerth, 2013/01/17)

<li> New: All vector classes should now have a in_local_range()
function indicating whether a given index is locally stored or not.
<br>
(Daniel Arndt, 2013/01/14)

<li> Fixed: The one-argument version of ConstraintMatrix::condense was
not prepared to deal with parallel vectors. This is now fixed.
<br>
(Daniel Arndt, Wolfgang Bangerth, 2013/1/9)

<li> New: Utilities::int_to_string can now deal with any 32-bit
integer, not just those with up to 6 digits.
<br>
(Wolfgang Bangerth, 2013/1/5)

<li> New: The PETScWrappers::MatrixBase::write_ascii() now takes a
(defaulted) argument allowing to select the PETSc Viewer style.
<br>
(Fahad Alrashed, 2013/1/2)

<li> New: The PETScWrappers::MPI::Vector::print() function overloads the
function of same name in the base class to ensure that the output
generated by a parallel vector makes sense.
<br>
(Fahad Alrashed, 2013/1/2)

<li> New: The PETScWrappers::VectorBase class now has a function
PETScWrappers::VectorBase::write_ascii() that allows writing the
vector's data to the default output stream.
<br>
(Fahad Alrashed, 2013/1/2)

<li> Fixed: PETScWrappers::SparseDirectMUMPS forgot to release
its memory upon destruction. This is now fixed.
<br>
(Alexander Grayver, 2012/12/12)

<li> Fixed: Using the copy constructor of FESystem led to trouble
down the road because some pointers were freed by the copy while
still in use by the original object. This is now fixed.
<br>
(Timo Heister, Wolfgang Bangerth, 2012/12/03)

<li> Improved: GridTools::make_periodicity_constraints substantially:
The low level interface now allows to specify a face orientation that
will be used when matching and constraining DoFs for periodic boundary
conditions. With that, the high level interface will work correctly on
(mainly 3d) grids that have cells in non standard orientation.
<br>
(Matthias Maier, 2012/12/01)

<li> Fixed: Fix GeometryInfo<2>::child_cell_on_face to respect face_flip
<br>
(Matthias Maier, 2012/12/01)

<li> New: There is now a version of DoFTools::make_zero_boundary_constraints()
that accepts a boundary indicator as argument.
<br>
(Wolfgang Bangerth, 2012/11/25)

<li> Fixed: The DoFTools::make_flux_sparsity_pattern() function
had a bug that triggered in 1d whenever there were neighboring
cells that differed in refinement level by more than one. This
is now fixed.
<br>
(Wolfgang Bangerth, 2012/11/20)

<li> Improved: The inner product, norm, and mean value computation
of deal.II's own Vector class are now parallelized by threads.
The parallelization does not change the order in which the additions
take place, ensuring exact reproducibility from one run to the next.
<br>
(Martin Kronbichler, 2012/11/18)

<li> New: The TrilinosWrappers::PreconditionBase class now has
a function TrilinosWrappers::PreconditionBase::Tvmult that
allows applying the transpose preconditioner.
<br>
(Guido Kanschat, 2012/11/04)

<li> New: The parallel::distributed::Triangulation::n_global_levels()
function returns the maximal refinement level over all involved
processors.
<br>
(Timo Heister, 2012/11/04)

<li> New: In addition to the regular subdomain ids (see
@ref GlossSubdomainId) there is now a second set of flags called
"level subdomain ids" that also assigns a subdomain to every cell
in a multigrid hierarchy.
<br>
(Timo Heister, 2012/11/04)

<li> Changed: GridOut::write_xfig has been improved in a number
of ways. In particular, one can now color cells based on a number
of different criteria that can be set in GridOutFlags::XFig.
<br>
(Guido Kanschat, 2012/11/04)

<li> The class Utilities::MPI::MPI_InitFinalize now also initializes
PETSc, when PETSc is installed.
<br>
(Timo Heister, 2012/11/02)

<li> Fixed: DoFTools::make_flux_sparsity_pattern wasn't prepared to
deal with adaptively refined meshes in 1d.
<br>
(Wolfgang Bangerth, 2012/10/30)

<li> New: Added PETScWrappers::PreconditionParaSails and
PETScWrappers::PreconditionNone. PETScWrappers::PreconditionParaSails
implements the interface to use the ParaSails sparse approximate
inverse preconditioner from the HYPRE suite. ParaSails supports
parallel distributed computations and can handle nonsymmetric
and also indefinite problems. PETScWrappers::PreconditionNone
implements non-preconditioning in PETSc which can be of use
together with the PETScWrappers::MatrixFree class.
<br>
(Martin Steigemann, 2012/10/26)

<li> New: The PETScWrappers::SparseDirectMUMPS class now allows to
exploit symmetry of the matrix, using the
PETScWrappers::SparseDirectMUMPS::set_symmetric_mode() function.
<br>
(Alexander Grayver, 2012/10/23)

<li> Fixed: Several static const member variables of the Accessor
classes were not properly instantiated. This only rarely created
trouble because they are typically only used as template arguments
and the compiler substituted them. However, one would get linker
errors when passing around a reference to them. This is now fixed.
<br>
(Wolfgang Bangerth, Guido Kanschat, 2012/10/11)

<li> Fixed: Handle lucky breakdowns in GMRES/FGMRES.
<br>
(B&auml;rbel Janssen, Timo Heister, 2012/10/09)

<li> Fixed: GridTools::find_cells_adjacent_to_vertex got into
trouble with anisotropically refined meshes. This is now fixed.
<br>
(Abner Salgado, Tobias Leicht, Wolfgang Bangerth, 2012/10/08)

<li> Fixed: FESystem can now deal with n_elements==0 for a block.
<br>
(Timo Heister, 2012/09/28)

<li> Fixed: ParameterHandler::print_parameters, when using
ParameterHandler::OutputStyle::LaTeX would always print a list
of parameters in each section as a latex itemized environment.
However, if there are none, we end up with an empty list which
latex does not like. This is now fixed.
<br>
(Wolfgang Bangerth, 2012/09/27)

<li> New: Added SparsityTools::distribute_sparsity_pattern() for
BlockCompressedSimpleSparsityPattern. This allows parallel computations
with distributed::Triangulation and PETScWrappers::MPI::BlockSparseMatrix.
<br>
(Timo Heister, 2012/09/25)

<li> New: Added BlockCompressedSimpleSparsityPattern::column_number().
<br>
(Timo Heister, 2012/09/25)

<li> New: There is now a function hp::FECollection::n_blocks() in analogy to
the existing function hp::FECollection::n_components().
<br>
(Wolfgang Bangerth, 2012/09/20)

<li> Changed: step-8 now outputs data in VTK format, rather than GMV.
GMV has long been dead.
<br>
(Wolfgang Bangerth, 2012/09/19)

<li> Fixed: One can compile deal.II with MPI support but run programs
that aren't intended to use parallel communications and that, in fact,
do not call <code>MPI_Init</code> at all. They are nevertheless supposed
to work but previously the TimerOutput would crash under these conditions.
This is now fixed.
<br>
(Timo Heister, Wolfgang Bangerth, 2012/09/18)

<li> Fixed: If you pipe content into the deallog object and there
is no end-line or flush after this content, and if a file stream
is associated to this object, and if that happens at the end of
the lifetime of the program, then the program would crash.
This is now fixed.
<br>
(Timo Heister, Wolfgang Bangerth, 2012/09/17)

<li> Fixed: The use of TableHandler::set_precision affected not only the
precision with which elements of a table were printed, but also the
precision carried by the output stream after writing the table was
finished. It thus affected the precision
with which later output was produced. This is now fixed.
<br>
(Timo Heister, 2012/09/16)

<li> Fixed: Output of super-columns in TableHandler::write_text()
was inconsistent. This is now fixed.
<br>
(Timo Heister, 2012/09/16)

<li> Changed: Due to incompatibilties with some hdf5 packages installed from
repositories we disable auto-detection of hdf5. Use --with-hdf if you need it.
<br>
(Timo Heister, 2012/09/14)

<li> New MeshWorker::LocalIntegrator and integration_loop() provide a
less confusing interface to MeshWorker loops.
<br>
(Guido Kanschat, 2012/09/13)

<li> New: TableHandler TextOutputFormat::simple_table_with_separate_column_description
that skips aligning the columns for increased performance.
<br>
(Timo Heister, 2012/09/10)

<li> Fixed: The Clang C++ compiler had some trouble dealing with obtaining
the return value of a Threads::Task object, due to a compiler bug in
dealing with friend declarations. This is now fixed.
<br>
(Wolfgang Bangerth, 2012/09/04)

<li> Fixed: When applying a ConstraintMatrix to a block matrix
where the last few rows are empty, we ran into an unrelated assertion.
This is now fixed.
<br>
(Jason Sheldon, Wolfgang Bangerth, 2012/09/04)
</ol>


*/
