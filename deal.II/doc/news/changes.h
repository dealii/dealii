/**
 * @page changes_after_6_3 Changes after Version 6.3

<p>
This is the list of changes made after the release of
deal.II version 6.3. It is subdivided into changes
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
  <li>None so far.
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
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
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li>None so far.
</ol>


<a name="lac"></a>
<h3>base</h3>

<ol>
  <li>
  <p>
  Fixed: deal.II release 6.3.0 did not compile with Trilinos versions 9.x and
  10.0. This is now fixed.
  <br>
  (Martin Kronbichler, WB 2010/06/28)
  </p>  
</ol>


<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
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
