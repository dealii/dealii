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
 @defgroup Iterators Iterators on mesh-like containers
 @{

deal.II has several classes which are understood conceptionally as
meshes. Apart from the obvious Triangulation, these are, for example,
DoFHandler and hp::DoFHandler. All of those define a set
of iterators, allowing the user to traverse the whole mesh, i.e. the
set of cells, faces, edges, etc that comprise the mesh, or portions of
it. These iterators are all in a sense derived from the TriaIterator
class.

Basically, the template signature of TriaIterator is
@code
  TriaIterator<Accessor>
@endcode

Conceptually, this type represents something like a pointer to an object
represented by the <code>Accessor</code> class.  Usually, you will not use the
actual class names spelled out directly, but employ one of the typedefs
provided by the container classes, such as <code>typename
Triangulation::cell_iterator</code>. Before going into this, let us
first discuss the concept of iterators, before delving into what the accessors
do.

As usual in C++, iterators, just as pointers, are incremented to the next
element using <tt>operator ++</tt>, and decremented to the previous element
using <tt>operator --</tt>. One can also jump <tt>n</tt> elements ahead using
the addition operator, <tt>it=it+n</tt>, and correspondingly to move a number
of elements back. In addition, and keeping with the tradition of the standard
template library, containers provide member functions <tt>begin()</tt> and
<tt>end()</tt> that provide the first element of a collection and a
one-past-the-end iterator, respectively. Since there are a number of different
iterators available, there is actually a whole family of such functions, such
as <tt>begin_active()</tt>, <tt>begin_face()</tt>, etc.

In terms of the concepts for iterators defined in the C++ standard, the
deal.II mesh iterators are bi-directional iterators: they can be incremented
and decremented, but an operation like <tt>it=it+n</tt> takes a compute time
proportional to <tt>n</tt>, since it is implemented as a sequence of
<tt>n</tt> individual unit increments. Note that this is in contrast to the
next more specialized iterator concept, random access iterators, for which
access to an arbitrary object requires only constant time, rather than linear.


@section IteratorsAndSets Iterators as pointers into sets of objects

As mentioned above, iterators in deal.II can be considered as iterating over
all the objects that constitute a mesh. (These objects are lines, quads, and
hexes, and are represented by the type of Accessor class given as template argument to the iterator.) This suggests to view a triangulation as a
collection of cells and other objects that are held together by a certain data
structure that links all these objects, in the same was as a linked list is
the data structure that connects objects in a linear fashion.

Triangulations in deal.II can indeed be considered in this way. In particular,
they use the computational notion of a forest of regular trees to store their
data. This can be understood as follows: Consider the cells of the coarse mesh
as roots; then, if one of these coarse mesh cells is refined, it will have
2<sup>dim</sup> children, which in turn can, but do not have to have
2<sup>dim</sup> children of their own, and so on. This means, that each cell
of the coarse mesh can be considered the root of a binary tree (in 1d), a
quadtree (in 2d), or an octree (in 3d). The collection of these trees
emanating from the cells of the coarse mesh then constitutes the forest that
completely describes the triangulation, including all of its active and
inactive cells. In particular, the active cells are those terminal nodes in
the tree that have no decendants, i.e. cells which are not further
refined. Correspondingly, inactive cells correspond to nodes in the tree with
descendents, i.e. cells that are further refined.

A triangulation contains forests for lines (each of which may have 2
children), quads (each with possibly four children), and hexes (each with no
or 8 children). Depending on the dimension, these objects are also termed
cells or faces.

Iterators loop over the elements of such forests. While the usual iterators
loop over all nodes of a forest, active iterators iterate over the
elements in the same order, but skip all non-active entries and therefore only
visit terminal nodes (i.e. active cells, faces, etc). There are many ways to
traverse the elements of a forest, for example breadth first or depth
first. Depending on the type of data structure used to store the forest, some
ways are more efficient than others. At present, the way iterators traverse
forests in deal.II is breadth first. I.e., iterators first visit all the
elements (cells, faces, etc) of the coarse mesh before moving on to all the
elements of the immediate level, i.e. the immediate children of the coarse
mesh objects; after this come the grandchildren of the coarse mesh, and so on.
However, it must be noted that programs should not rely on this particular
order of traversing a tree: this is considered an implementation detail that
can change between versions, even if we consider this an unlikely option at
the present time.



@section IteratorsDifferences Different kinds of iterators

Iterators have two properties: what they point to (i.e. the type of the
Accessor template argument), and the exact definition of the set they iterate
over. In general, iterators are always declared as
@code
  KindIterator<Accessor>
@endcode

Here, <tt>Kind</tt> determines what property an accessor needs to have to be
reached by this iterator (or omitted, for that matter). For example,
@code
  Iterator<Accessor>
@endcode
iterates over all objects of kind Accessor that make up the mesh (for example
all cells, whether they are further refined and have children, or not), whereas
@code
  ActiveIterator<Accessor>
@endcode
skips all objects that have children, i.e. objects that are not active.
Active iterators therefore operate on a subset of the objects
that normal iterators act on, namely those that possess the property that
they are active. Note that this is independent of the kind of object we
are operating on: all valid accessor classes have to provide the iterator
classes a method to find out whether they are active or not.

(For completeness, let us mention that there is a third kind of iterators: "raw
iterators" also traverse objects that are unused in the triangulation, but
allocated anyway for efficiency reasons. User code should never use raw
iterators, they are only for %internal purposes of the library.)

Whether an object is active can be considered a "predicate": a property that
is either true or false. Filtered iterators can be used to restrict the scope
of existing iterators even more. For instance, you could imagine to iterate
over the subset of those @ref GlossActive "active cells" having their user
flag set or belonging to a certain subdomain (both properties are either true
or false for a given object).

This is achieved by using an object of type FilteredIterator
&lt;BaseIterator&gt;, where BaseIterator usually is one of the
standard iterators discussed above.

The FilteredIterator gets an additional Predicate in its constructor and will
skip all objects where this Predicate evaluates to <tt>false</tt>. A
collection of predicates already implemented can be found in the namespace
IteratorFilters.


@subsection IteratorsLoops Iterating over objects

All iterators of the same kind and iterating over the
same kind of geometrical objects traverse the mesh in the same
order. Take this code example:
@code
  Triangulation<dim> tria;
  DoFHandler<dim>    dof1(tria);
  DoFHandler<dim>    dof2(tria);
  ...
  typename Trianguation<dim>::cell_iterator ti  = tria.begin();
  typename DoFHandler<dim>::cell_iterator   di1 = dof1.begin();
  typename DoFHandler<dim>::cell_iterator   di2 = dof2.begin();
  ...
  while (ti != tria.end())
  {
    // do something
    ++ti;
    ++di1;
    ++di2;
  }
@endcode

Here, all iterators will always point to the same mesh cell, even though
<tt>DoFHandler</tt> and <tt>Triangulation</tt> are very different classes,
and even if the DoFHandlers are handling different finite elements: they
all access cells in the same order, the difference is only in the Accessor.
As mentioned above, the order in which iterators traverse the forest of
objects is actually well-defined, but application programs should not
assume any such order, but rather consider this an implementation detail
of the library.

Corresponding to above example, the order in which iterators traverse active
objects is the same for all iterators in the following snippet, the difference to the previous example being that here we only consider active cells:
@code
  typename Trianguation<dim>::active_cell_iterator ti  = tria.begin_active();
  typename DoFHandler<dim>::active_cell_iterator   di1 = dof1.begin_active();
  typename DoFHandler<dim>::active_cell_iterator   di2 = dof2.begin_active();
  ...
  while (ti != tria.end())
  {
    // do something
    ++ti;
    ++di1;
    ++di2;
  }
@endcode



@section IteratorsAccessors Accessors

Iterators are like pointers: they can be incremented and decremented, but they
are really rather dumb. Their magic only lies in the fact that they point to
some useful object, in this case the Accessor. For pointers, they point to an
actual object that stores some data. On the other hand, the deal.II iterators,
when dereferenced, do not return a reference to an actual object, but return
an object that knows how to get at the data that represents cells. In general, this
object doesn't store itself where the vertices of a cell are or what its neighbors
are. However, it knows how to tease this sort of information from out of the
arrays and tables and lists that the Triangulation class sets up to describe a
mesh.

Accessing data that characterizes a cell is always done through the Accessor,
i.e. the expression <code>i-&gt;xxx()</code> grants access to <b>all</b>
attributes of this Accessor. Examples of properties you can query from an
iterator are
@code
  cell->vertex(1);
  line->child(0);
  hex->face(3);
  cell->at_boundary();
  face->boundary_indicator();
@endcode

Since dereferencing iterators yields accessor objects, these calls are to
member functions <code>Accesor::vertex()</code>,
<code>Accessor::child()</code> etc. These in turn figure out the relevant data
from the various data structures that store this data. How this is actually
done and what data structures are used is not really of concern to authors of
applications in deal.II. In particular, by hiding the actual data structures
we are able to store data in an efficient way, not necessarily in a way that
makes it easily accessible or understandable to application writers.



@section IteratorsTypedefs Kinds of accessors

Depending on what sort of data you want to access, there are different kinds
of accessor classes:

- The TriaAccessor class provides you with data that identifies the geometric
  properties of cells, faces, lines, quads, and hexes that make up a
  triangulation, as well as mother-child relationships.

- The CellAccessor class is derived from the TriaAccessor class for cases
  where an object has full dimension, i.e. is a cell rather than for example a
  line bounding a cell. In that case, additional information about the
  topological connection of a mesh is available from an accessor such as to
  request iterators pointing to neighbors of a cell.

- The DoFAccessor class lets you access information related to degrees
  of freedom associated with cells, faces, etc; it does so for both
  DoFHandler and hp::DoFHandler objects. Note that the DoFAccessor
  class is derived from either TriaAccessor or CellAccessor (depending
  on whether the DoFAccessor points to an object of full dimension or
  not) and so is able to provide a superset of information over its
  base classes. Additionally, the DoFAccessor class comes in two
  flavors, one accessing degrees of freedom on the level of a cell and
  the other accessing the active dofs of an active cell.

- The DoFCellAccessor class has the same purpose and relation to
  DoFCellAccessor as the CellAccessor has to TriaAccessor.

Except to look up member documentation, you will not usually have to deal with
the actual class names listed above. Rather, one uses the typedefs provided by
the container classes Triangulation, DoFHandler and hp::DoFHandler, as well
as the function that generate such objects:

<table border=1>
  <tr>
    <th>Container</th>
    <th>cell_iterator type</th>
    <th>function call</th>
  </tr>

  <tr>
    <th>Triangulation</th>
    <td>typename Triangulation::cell_iterator</td>
    <td>triangulation.begin()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>typename DoFHandler::cell_iterator</td>
    <td>dof_handler.begin()</td>
  </tr>

  <tr>
    <th>hp::DoFHandler</th>
    <td>typename hp::DoFHandler::cell_iterator</td>
    <td>hp_dof_handler.begin()</td>
  </tr>
</table>


<table border=1>
  <tr>
    <th>Container</th>
    <th>face_iterator type</th>
    <th>function call</th>
  </tr>

  <tr>
    <th>Triangulation</th>
    <td>typename Triangulation::face_iterator</td>
    <td>triangulation.begin_face()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>typename DoFHandler::face_iterator</td>
    <td>dof_handler.begin_face()</td>
  </tr>

  <tr>
    <th>hp::DoFHandler</th>
    <td>typename hp::DoFHandler::face_iterator</td>
    <td>hp_dof_handler.begin_face()</td>
  </tr>
</table>


Likewise, active iterators have the following properties:

<table border=1>
  <tr>
    <th>Container</th>
    <th>cell_iterator type</th>
    <th>function call</th>
  </tr>

  <tr>
    <th>Triangulation</th>
    <td>typename Triangulation::active_cell_iterator</td>
    <td>triangulation.begin_active()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>typename DoFHandler::active_cell_iterator</td>
    <td>dof_handler.begin_active()</td>
  </tr>

  <tr>
    <th>hp::DoFHandler</th>
    <td>typename hp::DoFHandler::active_cell_iterator</td>
    <td>hp_dof_handler.begin_active()</td>
  </tr>
</table>


<table border=1>
  <tr>
    <th>Container</th>
    <th>face_iterator type</th>
    <th>function call</th>
  </tr>

  <tr>
    <th>Triangulation</th>
    <td>typename Triangulation::active_face_iterator</td>
    <td>triangulation.begin_active_face()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>typename DoFHandler::active_face_iterator</td>
    <td>dof_handler.begin_active_face()</td>
  </tr>

  <tr>
    <th>hp::DoFHandler</th>
    <td>typename hp::DoFHandler::active_face_iterator</td>
    <td>hp_dof_handler.begin_active_face()</td>
  </tr>
</table>


In addition to these types and calls that act on cells and faces (logical
concepts that depend on the dimension: a cell is a quadrilateral in 2d, but
a hexahedron in 3d), there are corresponding types and calls like
<code>begin_active_quad()</code> or <code>end_quad()</code> that act on the
dimension independent geometric objects line, quad, and hex. These calls,
just as the ones above, exist in active and non-active forms.

The actual definition of all the typedefs local to the container classes are
stated in the

- internal::Triangulation::Iterators<1,spacedim>,
  internal::Triangulation::Iterators<2,spacedim>, and
  internal::Triangulation::Iterators<3,spacedim> classes for Triangulation
  iterators,

- internal::DoFHandler::Iterators<DH<1,spacedim> >,
  internal::DoFHandler::Iterators<DH<2,spacedim> >, and
  internal::DoFHandler::Iterators<DH<3,spacedim> > classes for DoFHandler
  and hp::DoFHandler iterators,


@section IteratorAccessorInternals Iterator and accessor internals

Iterators, being like pointers, act as if they pointed to an actual
object, but in reality all they do is to return an accessor when
dereferenced. The accessor object contains the state, i.e. it knows
which object it represents, by storing for example which Triangulation
it belongs to, and the level and index within this level of a cell. It
is therefore able to access the data that corresponds to the cell (or
face, or edge) it represents

There is a representation of past-the-end-pointers, denoted by special
values of the member variables <code>present_level</code> and <code>present_index</code> in
the TriaAccessor class: If <code>present_level</code> @> =0 and <code>present_index</code> @> =0,
then the object is valid; if
<code>present_level</code>==-1 and <code>present_index</code>==-1, then the iterator points
past the end; in all other cases, the iterator is considered invalid.
You can check this by calling the TriaAccessorBase::state() function.

Past-the-end iterators may also be used to compare an iterator with
the before-the-start value, when running backwards. There is no
distinction between the iterators pointing past the two ends of a
vector.

Cells are stored based on a hierarchical structure of levels, therefore the
above mentioned structure is useful. Faces however are not organized in
levels, and accessors for objects of lower dimensionality do not have a
<code>present_level</code> member variable.


@ingroup grid
*/

//@}


/**
 * @defgroup Accessors Accessor classes of the mesh iterators
 * @ingroup Iterators
 */


