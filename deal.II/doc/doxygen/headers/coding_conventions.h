//-------------------------------------------------------------------------
//    $Id: fe.h 12032 2006-01-15 00:41:50Z wolf $
//    Version: $Name$
//
//    Copyright (C) 1998, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @page CodingConventions Coding conventions used throughout deal.II
 
Throughout deal.II, we strive to keep our programming style and the kind
of interfaces we provide as consistent as possible. To this end, we
have adopted a set of coding conventions that we attempt to follow
wherever possible. When reading through them, it is important to
remember that they are not god given or better than any other set of
conventions; their purpose is merely to keep deal.II as uniform as
possible. Uniformity reduces the number of bugs we produce because we
can, for example, always assume that input arguments come before
output arguments of a function call. They also simplify reading code
because some things become clear already by looking at the style a
piece of code is written, without having to look up the exact
definition of something.

The following are the general rules we attempt to follow:
 
<ol>
<li> %Functions which return the number of something (number of cells,
  degrees of freedom, etc) should start with <code>n_*</code></li>

<li> %Function which set a bit or flag should start with <code>set_*</code>;
  functions which clear bits of flags should be named <code>clear_*</code></li>

<li> After each function, at least three empty lines are expected to
  enable better readability. One empty line occurs in functions to
  group blocks of code, two empty lines are not enough to distinguish
  visibly enough.</li>

<li> Whenever an integer variable can only assume nonnegative values,
  it has to be marked as unsigned.</li>

<li> Whenever an argument will not be changed, it should be marked
  const, even if it passed by value. This makes programs more readable
  and lets the compiler issue warnings if such a parameter variable is
  changed, which is often either involuntarily or poor style.</li>

<li> Whenever a function does not change any of the member variable of
  the embedding class/object, it should be marked as const.</li>

<li> %Function and variable names may not consist of only one or two
  letters, unless the variable is a pure counting index.</li>

<li> Use the geometry information in GeometryInfo<dim> to get the
  number of faces per cell, the number of children per cell, the
  child indices of the child cells adjacent to face 3, etc, rather
  than writing them into the directly as <code>2*dim</code>, <code>(1@<@<dim)</code> and
  <code>{0,3}</code>. This reduces the possibilities for errors and enhances
  readability of code. Unfortunately, the GeometryInfo mechanism
  was not invented right at the start of the program, so there are
  quite a lot of places where this rule is violated. Of you find
  such a place, fix it. We have set an amount of 1 cent per fixed
  place, payable by cheque when the deal.II project is considered
  finished by the author(s).</li>

<li> The layout of class declarations is the following: first the
  block of public functions, beginning with the constructors, then
  the destructors. If there are public member variables, these have
  to occur before the constructor. Public variables shall only be
  used if constant or unavoidable.</li>

  After the public members, the protected and finally the private
  members are to be listed. The order is as above: first variables
  then functions.</li>

  Exceptions shall be declared at the end of the public section
  before the non-public sections start.</li>

<li> If a function has both input and output parameters, usually the
  input parameters shall precede the output parameters, unless there
  are good reasons to change this order.</li>

<li> Exceptions are used for internal parameter checking and for
  consistency checks through the Assert macro. Exception handling
  like done by the C++ language (<code>try/throw/catch</code>) are used to
  handle run time errors (like I/O failures) which must be on
  in any case, not only in debug mode.</li>

<li> Classes and types generally are named using uppercase letters to denote
  word beginnings (e.g. TriaIterator) while functions and variables
  use lowercase letters and underscores to denote different words.
  The only exception are the iterator typedefs in Triangulation
  and DoFHandler (named cell_iterator, active_line_iterator, etc)
  to make the connection to the STL classes clear.</li>

<li> Each class has to have at least 200 pages of documentation ;-)</li>

</ol>
 
 */
