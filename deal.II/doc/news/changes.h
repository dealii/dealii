/**
 * @page changes_after_7_0 Changes after Version 7.0

<p>
This is the list of changes made after the release of
deal.II version 7.0.0.
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
  <li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
<li> Fixed: Boundary conditions in the step-23 tutorial program are now 
applied correctly. Matrix columns get eliminated with the used method
and introduce some contribution to the right hand side coming from
inhomogeneous boundary values. The old implementation did not reset the
matrix columns before applying new boundary values.<br>
(Martin Stoll, Martin Kronbichler, 2011/01/14)
</ol>

<ol>
<li> Extended: <code>base/tensor.h</code> has an additional collection of
contractions between three tensors (<i>ie</i>. <code>contract3</code>).
This can be useful for writing matrix/vector assembly in a more compact
form than before.<br>
(Toby D. Young, 2011/01/12)
</ol>

<a name="specific"></a>
<h3>Specific improvements</h3>


<ol>
  <li>
</ol>


*/
