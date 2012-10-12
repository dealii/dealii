/**
 * @page changes_after_7_2 Changes after Version 7.2

<p>
This is the list of changes made after the release of
deal.II version 7.2.0.
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

<li> Changed: the optional argument offset got removed from
DoFHandler and MGDoFHandler::distribute_dofs() because it was
never working correctly and it is not used.
<br>
(Timo Heister, 2012/09/03)
</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>


<ol>
<li> New: Added SparsityTools::distribute_sparsity_pattern() for
BlockCompressedSimpleSparsityPattern. This allows parallel computations
with distributed::Triangulation and PETScWrappers::MPI::BlockSparseMatrix.
<br>
(Timo Heister, 2012/09/25)

<li> Simplifications of the internal structures of Triangulation and
DoFHandler, in particular removal of specializations.
<br>
(Guido Kanschat, 2012/09/13)
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
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
