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
 * @page changes_between_3_2_0_and_3_2_1 Changes between Version 3.2.0 and 3.2.1

<p>
This is the list of changes made between the deal.II releases listed above.
All entries are signed with the names of the author.
</p>


<h2>Changes between versions 3.2.0 and 3.2.1</h2>

<ol>
  <li> <p>
       Fixed: In the <code>ParameterHandler</code>
       class, we leaked 8 or 12 bytes of memory per declared
       parameter. This is now fixed.
       <br>
       (WB 2001/11/28)
       </p>

  <li> <p>
       Fixed: he <code>DoFHandler</code> class had a
       memory leak. This is now fixed.
       <br>
       (WB 2001/11/28)
       </p>

  <li> <p>
       Fixed: The <code
       class="member">FilteredMatrix::allocate_tmp_vector</code>
       function had a bug with block vectors.
       <br>
       (WB 2001/11/22)
       </p>

  <li> <p> New: example program step-12. Discontinuous Galerkin
       discretization.
       <br>
       (RH 2001/11/21)
       </p>

  <li> <p> New: The new <code
       class="member">CellAccessor<dim>::neighbor_of_coarser_neighbor</code>
       function returns where to find the present cell from a coarser
       neighbor.
       <br>
       (RH 2001/11/21)
       </p>

  <li> <p>
       Fixed: when checking for convergence in linear solvers in
       <code>SolverControl::check</code>, we first
       checked whether the maximal iteration count was reached, and
       only then whether the target residual was achieved. In cases,
       where the target residual was only reached in the very last
       allowed iteration, this led to a failure notice of the linear
       solver, rather than to a success message. This is now fixed.
       <br>
       (WB 2001/11/19)
       </p>

  <li> <p>
       Fixed: an error in the definition of the <code
       class="member">SolverMinRes::solve</code> function prevented
       its compilation under some circumstances.
       <br>
       (WB 2001/11/14)
       </p>

  <li> <p>
       Fixed: upon breakdown, the <code
       class="class">SolverBicgstab</code> forgot to increment the
       iteration counter for the breakdown cycle. This is now fixed.
       <br>
       (WB 2001/11/14)
       </p>

  <li> <p>
       Fixed: the <code
       class="member">DoFTools::make_flux_sparsity_pattern</code> was
       implemented in 1d, but not instantiated. This is now fixed.
       <br>
       (WB 2001/11/14)
       </p>

  <li> <p>
       Fixed: the <code
       class="member">GridGenerator::hyper_rectangle</code> function
       was broken in 3d.
       <br>
       (WB 2001/10/26)
       </p>

  <li> <p>
       Fixed: class <code>SparsityPattern</code> can
       handle rows of length zero now. For quadratic matrices, these
       rows are changed to length one, since a diagonal element must
       be stored.
       <br>
       (GK 2001/10/11)
       </p>

  <li> <p>
       Fixed: The <code>DataOutBase::EpsFlags</code>
       class forgot to declare the reverse grey scale function as one
       possible input for the color function for the
       <code>ParameterHandler</code> class. This is now
       possible.
       <br>
       (WB 2001/10/10)
       </p>

  <li> <p>
       Fixed: the iterator category template base class of grid
       iterators was incorrectly set.
       <br>
       (WB 2001/09/28)
       </p>
</ol>


*/
