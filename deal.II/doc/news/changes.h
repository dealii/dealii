/**
 * @page changes_after_6.1 Changes after Version 6.1

<p>
This is the list of changes made after the release of 
deal.II version 6.1. It is subdivided into changes
made to the three sub-libraries <a href="#base">base</a>, 
<a href="#lac">lac</a>, and <a href="#deal.II">deal.II</a>, as well as
changes to the <a href="#general">general infrastructure,
documentation, etc</a>.
<br>
Furthermore, all changes which have been performed en route to
extending deal.II from isotropic to anisotropic refinement are listed
in this <a href="anisotropic_changes.html">list</a>.
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
  <li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li>
  <p>
  Fixed: In the new @ref step_33 "step-33" tutorial program there was
  a place where we incorrectly passed the diameter of a cell, rather
  than a face, to a function. This is now fixed.
  <br>
  (Chih-Che Chueh, WB 2008/06/28)
  </p>
  </li>

  <li>
  <p>
  Fixed: The new @ref step_33 "step-33" tutorial program had a place where
  we didn't release some memory. This is now fixed.
  <br>
  (Chih-Che Chueh, WB 2008/06/02)
  </p>
  </li>

  <li>
  <p>
  Fixed: A missing include file prevented the <code>./configure</code> script
  from detecting the presence of the demangler with recent versions of the
  gcc compiler. The result is that backtraces after failed assertions only
  show the mangles function names, not their plain text equivalent. This is
  now fixed.
  <br>
  (WB 2008/05/27)
  </p>
  </li>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li>
  <p>
  Fixed: When using a higher order mapping of degree at least 3 (i.e.
  the MappingQ class) on meshes that have cells with non-standard
  face orientation and that touch the boundary of the domain, then
  some interpolation points were wrongly computed.
  <br>
  (Tobias Leicht, Timo Heister, WB 2008/06/10)
  </p>
  </li>

  <li>
  <p>
  Fixed: The ParameterHandler::get_integer function now throws an
  exception when called for non-integer parameters. This exception has
  been commented out some time ago but is now reincluded.
  <br>
  (RH 2008/06/11)
  </p>
  </li>
</ol>


*/
