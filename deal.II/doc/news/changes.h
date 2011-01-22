/**
 * @page changes_after_7_0 Changes after Version 7.0

<p>
This is the list of changes made after the release of
deal.II version 7.0.0.
All entries are signed with the names of the author.
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
  <li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
<li> Fixed: When using Trilinos and using the Intel C++ compiler,
we accidentally used invalid compiler flags that led to a warning
every time we compiled a file..
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: At the bottom of the page of tutorial programs we show a "plain"
version of the tutorial program. However, the script that generates this plain
version was broken and sometimes truncated the file. This
should be fixed now.
<br>
(Wolfgang Bangerth, 2011/01/18)

<li> Extended: Several missing instantiations of functions for triangulations and DoF handlers embedded in higher dimensional space have been added.
<br>
(Wolfgang Bangerth, 2011/01/15)
</ol>



<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
<li> Changed: The TrilinosWrappers::MPI::BlockVector::compress function now takes an
argument (with a default value) in exactly the same way as the
TrilinosWrappers::MPI::Vector::compress function already did.
<br>
(Wolfgang Bangerth, 2011/01/21)
</li>

<li> Fixed: When calling DataOut::build_patches with a mapping, requesting more
than one subdivision, and when <code>dim@<spacedim</code>, then some cells
were not properly mapped. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/01/18)
</li>

<li> New: Restructured the internals of PETScWrappers::Precondition* to allow a
PETSc PC object to exist without a solver. New: use Precondition*::vmult() to
apply the preconditioner once. Preconditioners now have a default constructor
and an initialize() function and are no longer initialized in the solver call,
but in the constructor or initialize().
<br>
(Timo Heister, 2011/01/17)
</li>

<li> Fixed: Boundary conditions in the step-23 tutorial program are now
applied correctly. Matrix columns get eliminated with the used method
and introduce some contribution to the right hand side coming from
inhomogeneous boundary values. The old implementation did not reset the
matrix columns before applying new boundary values.<br>
(Martin Stoll, Martin Kronbichler, 2011/01/14)
</li>

<li> Extended: <code>base/tensor.h</code> has an additional collection of
contractions between three tensors (<i>ie</i>. <code>contract3</code>).
This can be useful for writing matrix/vector assembly in a more compact
form than before.<br>
(Toby D. Young, 2011/01/12)
</ol>


*/
