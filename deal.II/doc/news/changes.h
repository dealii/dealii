/**
 * @page changes_after_6.2 Changes after Version 6.2

<p>
This is the list of changes made after the release of 
deal.II version 6.2. It is subdivided into changes
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
  <li>
  <p> 
  </p>
  </li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>
  <li>
  <p>
  </p>
  </li>
</ol>



<a name="base"></a>
<h3>base</h3>

<ol>
  <li>
  <p>
  </p>
  </li>
</ol>



<a name="lac"></a>
<h3>lac</h3>

<ol>
  <li>
  <p>
  </p>
  </li>
</ol>



<a name="deal.II"></a>
<h3>deal.II</h3>

<ol>
  <li>
  <p>
  New: The DoFTools::make_sparsity_pattern functions have acquired a
  new paramater <code>subdomain_id</code>. If a value other than the
  default value is passed for it, the function only assembles the
  sparsity pattern on those cells that have the given subdomain id.
  This is useful, for example, in conjunction with the
  TrilinosWrappers::SparsityPattern class that can hold a sparsity
  pattern distributed across several MPI processes; in that case, it is
  not necessary that each process builds the entire sparsity pattern.
  <br>
  (WB 2009/04/29)
  </p>
  </li>
</ol>


*/
