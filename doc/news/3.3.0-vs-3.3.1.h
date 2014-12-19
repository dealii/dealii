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
 * @page changes_between_3_3_0_and_3_3_1 Changes between Version 3.3.0 and 3.3.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>


<ol>
  <li> <p> 
       Fixed: In 3d, the function <code
       class="member">DoFTools::make_hanging_node_constraints</code> 
       contained an assertion that failed erroneously for finite
       elements that do not have degrees of freedom on vertices. This
       is now fixed.
       <br> 
       (WB 2002/02/21)
       </p>

  <li> <p>
       Fixed: <code>TriaAccessor<3,3>::measure</code>
       sometimes computed a negative value. This is now fixed.
       <br> 
       (WB 2002/02/21)
       </p>
</ol>

*/
