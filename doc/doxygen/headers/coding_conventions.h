// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


/**
 * @page CodingConventions Coding conventions used throughout deal.II

Throughout deal.II, we strive to keep our programming style and the kind of
interfaces we provide as consistent as possible. To this end, we have adopted
a set of coding conventions that we attempt to follow wherever possible. They
have two parts: style issues, and something we call "defensive programming",
the latter being an attempt to let our code help us find bugs.  When reading
through them, it is important to remember that styles are not god given or
better than any other set of conventions; their purpose is merely to keep
deal.II as uniform as possible. Uniformity reduces the number of bugs we
produce because we can, for example, always assume that input arguments come
before output arguments of a function call. They also simplify reading code
because some things become clear already by looking at the style a piece of
code is written, without having to look up the exact definition of something.

<h3>Notes on deal.II indentation</h3>

<p>deal.II uses <code>clang-format</code> 6.0 to normalize indentation. A
style file is provided at
@code
  \${SOURCE_DIR}/.clang-format
@endcode

<p>Before a commit, you should run
@code
  clang-format -i <file>
@endcode
on each of your files. This will make sure indentation is conforming to the
style guidelines outlined in this page.

This is cumbersome. Consequently, and more easily, you can just run
@code
  make indent
@endcode
in whatever directory you set up the library to be compiled in, to indent all
source files that have been changed recently. If you want to make sure that
the indenting is correct for all your commits, you might want to set up a
pre-commit hook. One way to do so, is to copy
<code>\${SOURCE_DIR}/contrib/git-hooks/pre-commit</code> to
<code>\${SOURCE_DIR}/.git/hooks/pre-commit</code> and make sure it is
executable.

If the system you are working on has more than one version of
<code>clang-format</code> installed (or if it is not in the path)
you should replace the above <code>make indent</code> command with
@code
  make DEAL_II_CLANG_FORMAT=/path/to/clang-6.0/clang-format indent
@endcode
to point to the correct executable.
</p>

<h3>Style issues</h3>

<ol>
<li> %Functions which return the number of something (number of cells,
  degrees of freedom, etc) should start with <code>n_*</code>. Example:
  SparsityPatternBase::n_nonzero_elements().</li>

<li> %Functions which set a bit or flag should start with <code>set_*</code>;
  functions which clear bits or flags should be named <code>clear_*</code>.
  Example: CellAccessor::set_refine_flag().</li>

<li> Traditional logical operators should be used instead of their English
  equivalents (i.e., use <code>&&</code>, <code>||</code>, and <code>!</code>
  instead of <code>and</code>, <code>or</code>, and <code>not</code>).

<li> In the implementation files, after each function, three empty lines are
  expected to enable better readability. One empty line occurs in functions to
  group blocks of code, since two empty lines are not enough to visibly
  distinguish sufficiently that the code belongs to two different
  functions.</li>

<li> Whenever an integer variable can only assume nonnegative values,
  it is marked as unsigned. The same applies to functions that can only
  return positive or zero values. Example: Triangulation::n_active_cells().</li>

<li> Whenever an argument to a function will not be changed, it should be marked
  const, even if passed by value. Generally, we mark input parameters as
  const. This aids as an additional documentation tool to clarify the
  intent of a parameter (input, output, or both)
  and lets the compiler issue warnings if such a parameter is
  changed, which is often either involuntarily or poor style.</li>

<li> Whenever a function does not change any of the member variable of
  the embedding class/object, it should be marked as const.</li>

<li> %Function and variable names may not consist of only one or two
  letters, unless the variable is a pure counting index.</li>

<li> Type aliases (<code>using</code>-declarations) are preferred to
  <code>typedef</code>-declarations.</li>

<li> Use the geometry information in GeometryInfo to get the
  number of faces per cell, the number of children per cell, the
  child indices of the child cells adjacent to face 3, etc, rather
  than writing them into the code directly as <code>2*dim</code>,
  <code>(1@<@<dim)</code> and
  <code>{0,3}</code>. This reduces the possibilities for errors and enhances
  readability of code.</li>

<li> The layout of class declarations is the following: first the
  block of public functions, beginning with the constructors, then
  the destructors. If there are public member variables, these have
  to occur before the constructor. Public variables shall only be
  used if constant (in particular if they are static and constant)
  or unavoidable.
  <br>
  After the public members, the protected and finally the private
  members are to be listed. The order is as above: first variables
  then functions.
  <br>
  Exceptions shall be declared at the end of the public section
  before the non-public sections start.
  <br>
  We do not use the C++11-style class member initialization for member variables
  that are neither <code>static const</code> nor <code>static constexpr</code>;
  i.e., instead of
@code
  class Foo
  {
    int a = 42;
    int *b = nullptr;
  };
@endcode
  write
@code
  class Foo
  {
    Foo();

    int a;
    int *b;
  };



  inline Foo::Foo()
  : a(42)
  , b(nullptr)
  {}
@endcode
  </li>

<li> If a function has both input and output parameters, usually the
  input parameters shall precede the output parameters, unless there
  are good reasons to change this order. (The most common reason is trailing
  input parameters with default values.) </li>

<li> Exceptions are used for %internal parameter checking and for
  consistency checks through the Assert macro. Exception handling
  like done by the C++ language (<code>try/throw/catch</code>, and using the
  AssertThrow macro) are used to
  handle run time errors (like I/O failures) which must be on
  in any case, not only in debug mode.</li>

<li> Sometimes it makes sense to implement a class by using several
  non-member functions that are not part of the public interface and are only
  meant to be called in the current source file. Such free functions should be
  put in an internal namespace structured in the following way:
  @code
  namespace internal
  {
    namespace ClassNameImplementation
    {
      // free functions go here
    }
  }
  @endcode
  where <code>ClassName</code> is the name of the calling class.

<li> Classes, namespaces and types generally are named using uppercase letters
  to denote word beginnings (e.g. TriaIterator) &mdash; sometimes called
  <a href="http://en.wikipedia.org/wiki/Camel_case"><i>camel
  case</i></a> &mdash; while functions and variables
  use lowercase letters and underscores to separate words.
  The only exception are the iterator alias in Triangulation
  and DoFHandler (named cell_iterator, active_line_iterator, etc)
  to make the connection to the standard library container classes clear.</li>

<li> For classes with multiple template arguments, the dimension is usually
  put before the data type specifier, i.e., we use Point<dim,number> and not
  Point<number,dim>.

<li> There are several places in deal.II where we use forward declarations in
  header files. The reason for this is that we can, hopefully, improve
  compilation speeds by not using headers when we just need to mark a certain
  type as an argument to a function. The convention used in deal.II is that, if
  all we need is a type name, then the type may be forward declared in the
  header where we need it; if a function (or member function) can return a value
  then a declaration of that value's type should be available (by including the
  necessary header). For example, <code>deal.II/dofs/dof_handler.h</code>
  includes <code>deal.II/dofs/dof_accessor.h</code> so that one can write
  something like <code>dof_handler.begin_active()->is_active()</code> without
  explicitly including the header declaring the type of the object returned by
  <code>begin_active()</code>.

<li> Each class has to have at least 200 pages of documentation ;-)</li>

</ol>


<h3>Instantiation of templated functions/classes</h3>

<p>The majority of classes and functions in deal.II are templated. This brings a
question of how and where such objects are instantiated, if at all. Throughout
deal.II we adopt the following convention:</p>

<ol>

<li> If we can enumerate all possible template arguments (e.g., the dimension
can only be 1, 2, or 3), then a function template goes into the <code>.cc</code>
file and we explicitly instantiate all possibilities. Users will not have any
need to ever see these function templates because they will not want to
instantiate these functions for any other template arguments anyway. </li>

<li> If we can not enumerate all possible template arguments (e.g., vector
types -- because users might want to define their own vector kinds) but at
least know a few common usage cases, then the function is put into a
<code>.templates.h</code> file. We \#include it into the <code>.cc</code> file
and instantiate the functions for all of the common arguments. For almost all
users, this will be just fine -- they only use the (vector, matrix, ...) types
we already instantiate, and for them the <code>.templates.h</code> file will not
be of any interest. It will also not slow down their compilations because
nothing they see will \#include the <code>.templates.h</code> file. But users who
define their own (vector, matrix, ...) types can instantiate the template
functions with their own user-defined types by including the
<code>.templates.h</code> files.

<li> Finally, if we can not assume in advance which values template arguments
will take (e.g., any class derived from Subscriptor can be used as an argument),
the definitions of functions are provided at the bottom of the header
file with declarations. The definitions should be guarded with <code>\#ifndef
DOXYGEN ... \#endif</code> to prevent Doxygen from picking them up.</li>

</ol>

<p> For the first two cases, instantiation instructions are defined in
<code>.inst.in</code> files. They are processed by a binary called
expand_instantiations (built from
<code>cmake/scripts/expand_instantiations.cc</code>) and the parameters are
defined dynamically through cmake depending on your configuration (see
<code>cmake/config/template-arguments.in</code> in your build directory).
It is those <code>.inst</code> files that are eventually included from the
corresponding <code>.cc</code> files. </p>


<h3>Defensive programming</h3>

<p> Defensive programming is a term that we use frequently when we talk about
writing code while in the mindset that errors will happen. Here, errors can
come in two ways: first, I can make a mistake myself while writing a
functions; and secondly, someone else can make a mistake while calling my
function. In either case, I would like to write my code in such a way that
errors are (i) as unlikely as possible, (ii) that the compiler can already
find some of the mistakes, and (iii) that the remaining mistakes are
relatively easy to find, for example because the program aborts. Defensive
programming is then a set of strategies that make these goals more likely.
</p>

<p>
Over time, we have learned a number of techniques to this end, some of which
we list here:
<ol>
<li> <i>Assert preconditions on parameters:</i> People call functions with wrong
  or nonsensical parameters, all the time. As the prototypical example,
  consider a trivial implementation of vector addition:
  @code
    Vector &
    operator+=(Vector       &lhs,
               const Vector &rhs)
    {
      for (unsigned int i=0; i<lhs.size(); ++i)
        lhs(i) += rhs(i);
      return lhs;
    }
  @endcode
  While correct, this function will get into trouble if the two vectors
  do not have the same size. You think it is silly to call this function
  with vectors of different size? Yes, of course it is. But it happens
  all the time: people forget to reinitialize a vector, or it is reset in
  a different function, etc. It happens. So if you are in such an unlucky
  case, it can take a long time to figure out what's going on because
  you are likely to just read uninitialized memory, or maybe you are
  writing to memory the <code>lhs</code> vector doesn't actually own.
  Neither is going to lead to immediate termination of the program,
  but you'll probably get random errors at a later time. It would be
  much easier if the program just stopped here right away. The following
  implementation will do exactly this:
  @code
    Vector &
    operator+=(Vector       &lhs,
               const Vector &rhs)
    {
      Assert (lhs.size() == rhs.size(),
              ExcDimensionMismatch(lhs.size(), rhs.size());
      for (unsigned int i=0; i<lhs.size(); ++i)
        lhs(i) += rhs(i);
      return lhs;
    }
  @endcode
  The <code>Assert</code> macro ensures that the condition is true
  at run time, and otherwise prints a string containing information
  encoded by the second argument and aborts the program. This way,
  when you write a new program that happens to call this function,
  you will learn of your error right away and have the opportunity
  to fix it without ever having to seriously debug anything.
  <p>
  As a general guideline, whenever you implement a new function,
  think about the <i>preconditions</i> on parameter, i.e. what does the
  function expect to be true about each of them, or their combination.
  Then write assertions for all of these preconditions. This may be
  half a dozen assertions in some cases but remember that each assertion
  is a potential bug already found through trivial means.
  <p>
  In a final note, let us remark that assertions are of course expensive:
  they may make a program 3 or 5 times slower when you link it against
  the debug version of the library. But if you consider your <i>overall</i>
  development time, the ability to find bugs quickly probably far outweighs
  the time you spend waiting for your program to finish. Furthermore,
  calls to the Assert macro are removed from the program in optimized mode
  (which you presumably only use once you know that everything runs just
  fine in debug mode. The optimized libraries are faster by a factor of
  3-5 than the debug libraries, at the price that it's much harder to find
  bugs.
  </li>

<li> <i>Assert postconditions:</i> If a function computes something
  non-trivial there may be a bug in the code. To find these, use
  postconditions: just like you have certain knowledge about useful values
  for input parameters, you have knowledge about what you expect possible
  return values to be. For example, a function that computes the norm of
  a vector would expect the norm to be positive. You can write this as
  follows:
  @code
    double norm(const Vector &v)
    {
      double s = 0;
      for (unsigned int i=0; i<v.size(); ++i)
        s += v(i) * v(i);

      Assert (s >= 0, ExcInternalError());
      return std::sqrt(s);
    }
  @endcode
  This function is too simple to really justify this assertion, but imagine
  the computation to be lengthier and you can see how the assertion helps
  you ensure (or <i>hedge</i>) yourself against mistakes. Note that one
  could argue that the assertion should be removed once we've run the program
  a number of times and found that the condition never triggers. But it's
  better to leave it right where it is: it encodes for the future (and for
  readers) knowledge you have about the function; if someone comes along
  and replaced the implementation of the function by a more efficient
  algorithm, the assertion can help make sure that the function continues
  to do what it is supposed to do.
  </li>

<li> <i>Assert internal states:</i> In a similar vein, if you have a
  complex algorithm, use assertions to ensure that your mental model
  of what is going on matches what is indeed true. For example, assume
  you are writing a function that ensures that mesh sizes do not change
  too much locally. You may end up with a code of the following kind:
  @code
    for (const auto &cell = triangulation.active_cell_iterators())
      for (unsigned int face=0; ...)
        {
          if (something)
            { ... }
          else
            {
              // we have a cell whose neighbor must
              // be at the boundary if we got here
            }
        }
  @endcode
  The conditions that got us into the else-branch may be
  complicated, and while it may be true that we believed that the
  only possibility we got here is that the neighbor is at the boundary,
  there may have been a bug in our implementation. There may also have been
  a bug in our thinking, or someone changes the code way above in the same
  function and forgets about the issue here, or a change at a completely
  different location in the library makes the assumption untenable. In
  all of these cases, the explicit statement of our assertion makes sure
  that these problems are easily found.
  </li>

<li> <i>Initialize variables at the point of their declaration if they
  live on the stack:</i>
  Traditional C required that variables are declared at the beginning of
  the function even if they are only used further below. This leads to
  code like this that we may imagine in a 1d code:
  @code
    template <int dim>
    void foo ()
    {
      Point<dim> cell_center;
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
          ...
        }
      ...
    }
  @endcode
  The problem is that if the code between the declaration and initialization
  is long and complicated, you can't look up on one page what the type of
  a variable is and what it's value may be. In fact, it may not even be
  quite clear that the variable is used initialized at all, or whether it
  is accidentally left uninitialized.
  <p>
  A better way to do this would be as follows:
  @code
    template <int dim>
    void foo ()
    {
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          Point<dim> cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
          ...
        }
      ...
    }
  @endcode
  This makes it much clearer what the type of the variable is
  and that it is in fact only ever used when initialized. Furthermore,
  if someone wants to read the code to see what the variable is in fact
  doing, declaring and initializing it in the innermost possible scope
  makes this task easier: we don't have to look upwards for it beyond
  the declaration, and we don't have to look downward beyond the end
  of the current scope since this is where the variable dies.
  <p>
  As a final note, it is clear that you can only do this sort of stuff
  for variables that completely live on the stack without allocating memory
  on the heap. Within deal.II, this is only true for builtin types like
  <code>int, double, char</code>, etc, as well as the Point and Tensor
  classes. Everything else has something like a <code>std::vector</code>
  as a member variable, which requires memory allocation &mdash; you don't
  want to declare these inside loops, at least not if the loop is
  traversed frequently.
  </li>

<li> <i>Make variables const:</i> To pick up on the example above, note
  that in most cases we will never change the variable so initialized
  any more. In other words, if this is the case, we may as well write
  things as follows:
  @code
    template <int dim>
    void foo ()
    {
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          const Point<dim> cell_center = (cell->vertex(0) +
                                          cell->vertex(1)) / 2;
          ...
        }
      ...
    }
  @endcode
  By marking the variable as constant we make sure that we don't accidentally
  change it. For example, the compiler could catch code like this:
  @code
        if (cell_center[0] = 0)
          ...
  @endcode
  This was most likely meant to be a <code>==</code> rather than an
  assignment. By marking the variable as const, the compiler would have
  told us about this bug. Maybe equally importantly, human readers of the
  code need not look further down whether the value of the variable may
  actually be changed somewhere between declaration and use &mdash; it
  can't be if it is marked as const.
  </li>

<li> <i>Make input arguments of functions const:</i> The same essentially
  holds true as well as for function arguments: If you have no intention
  of changing a variable (which is typically the case for input arguments),
  then mark it as constant. For example, the following function should take
  its argument as a constant value:
  @code
     template <int dim>
     typename Triangulation<dim>::cell_iterator
     CellAccessor<dim>::child(const unsigned int child_no)
     {
       ...
       return something;
     }
  @endcode
  Here, the user calls <code>cell-@>child(3)</code>, for example. There really
  is no reason why the function would ever want to change the value of
  the <code>child_no</code> argument &mdash; so mark it as constant:
  this both helps the reader of the code understand that this is an
  input argument of the function for which we need not search below whether
  it is ever changed, and it helps the compiler help us finding bugs if
  we ever accidentally change the value.
</ol>

 */
