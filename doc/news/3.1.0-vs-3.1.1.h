// ---------------------------------------------------------------------
//
// Copyright (C) 2013, 2014 by the deal.II authors
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
 * @page changes_between_3_1_0_and_3_1_1 Changes between Version 3.1.0 and 3.1.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>

<ol>
  <li> <p>
       Fixed: Templatized classes which had a default template
       argument that contains colons (such as nested types) did
       not show up in the forward declaration files, and were also
       missing from the class index. This is fixed now.
       <br>
       (WB 2001/05/02)
       </p>

  <li> <p>
       Fixed: The vertex numbers written by the <code
       class="class">GridOut</code>::<code
       class="member">write_ucd_faces</code> function are now also
       1-based. This, to be consistent with the vertex numbers given
       by the <code>GridOut</code>::<code
       class="member">write_ucd</code> function.
       <br>
       (RH 2001/04/20)
       </p>

  <li> <p>
       Fixed: the <code
       class="member">DoFRenumbering::Cuthill_McKee</code> function
       did not work correctly when giving the <code
       class="member">reversed_numbering</code> flag (off-by-one
       indexing). This is now fixed.
       <br>
       (<a href="mailto:or@winfos.com">Oliver Rheinbach</a> 2001/04/12)
       </p>
       
  <li> <p>
       Fixed: When using Neuman boundary functions in the 
       <code>KellyErrorEstimator</code> class, it was
       assumed that the function object had <code
       class="member">Function::vector_value</code> overloaded, even
       in the scalar case. We now use <code
       class="member">Function::value</code> instead.
       <br>
       (WB 2001/04/09)
       </p>
       
  <li> <p>
       New/Fixed: Now there exists a new <code
       class="class">Triangulation</code>::<code
       class="member">ExcMultiplySetLineInfoOfLine</code> exception,
       that is thrown if the <code>SubCellData</code>
       that is given to <code
       class="class">Triangulation</code>::<code
       class="member">create_triangulation</code>, multiply includes
       the line info of a specific line. Before the fix the wrong
       <code>ExcInteriorLineCantBeBoundary</code>
       exception was thrown.
       <br>
       (RH 2001/04/03)
       </p>
       
  <li> <p>
       Fixed: Missing <code>ucd</code> is now added to the list of
       supported output formats returned by <code
       class="class">GridOut</code>::<code
       class="member">get_output_format_names</code>.
       <br>
       (RH 2001/04/03)
       </p>

  <li> <p>
       Fixed: the program that generated HTML files from CVS back logs
       was run unintentionally. It failed because the distribution
       does not have the necessary CVS files. The program is no more
       run now.
       <br>
       (WB 2001/03/23)
       </p>

  <li> <p>
       Fixed: the program that generated HTML from the example
       programs was broken on some platforms for reasons beyond our
       knowledge. This is now fixed.
       <br>
       (Roger Young, WB 2001/03/22)
       </p>

  <li> <p> 
       Removed: The explicite instantiations of <code
       class="class">SparseMatrix&lt;long double&gt;</code> are removed as a
       prerelease of gcc3.0 fails to compile it. Now the user of <code
       class="class">SparseMatrix&lt;long double&gt;</code> needs to include
       <code>lac/include/lac/sparse_matrix.templates.h</code> into his
       source file and to use a different compiler, e.g. gcc 2.95.2 or
       a future version of gcc3.0 (that will then hopefully be fixed).
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       Changed: We now unconditionally include
       <code>deal.II/include/grid/tria_accessor.templates.h</code>
       (which contains some inline functions for triangulation
       accessor classes) into 
       <code>deal.II/include/grid/tria_accessor.h</code> to work
       around a problem with gcc3.0 which does not place a copy of
       these functions into the library. Previously we only included
       the file in optimized mode.
       <br>
       (RH 2001/03/14)
       </p>

  <li> <p>
       New: The class <code>GridReordering::Cell</code> has
       now a copy constructor to work around a bug in a gcc3.0
       snapshot.
       <br>
       (RH, WB 2001/03/14)
       </p>

  <li> <p>
       Fixed: a missing <code>const</code> prevented the <code
       class="class">PreconditionSelector</code> class to be
       compiled.
       <br>
       (Roger Young, John Burnell 2001/03/08)
       </p>

  <li> <p>
       Updated: the scripts <code>config.sub</code> and
       <code>config.guess</code> were updated to their newest
       versions, since <code>./configure</code> did not run properly
       on all systems. They are needed to determine the host system
       type correctly.
       <br>
       (WB 2001/03/06)
       </p>

  <li> <p>
       Fix: in the triangulation, the <code
       class="member">straight_boundary</code> variable, which is a
       reference, was assigned the address of a temporary object. It
       is unclear how this could have worked for three years, but it
       apparently did...
       <br>
       (WB 2001/02/26)
       </p>

  <li> <p>
       Fix: the <code
       class="member">DoFTools::compute_intergrid_constraints</code>
       function took memory quadratic in the number of degrees of
       freedom. This is now reduced to linear behaviour, with a
       constant that depends on the number of levels by which the two
       grids differ.
       <br>
       (WB 2001/02/26)
       </p>
</ol>


*/
