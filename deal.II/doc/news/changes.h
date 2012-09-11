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
<li> Nothing so far.
</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
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
