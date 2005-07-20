/**
 @defgroup Iterators Iterators on mesh-like containers
 @{

deal.II has several classes which are understood conceptionally as
meshes. Apart from the obvious Triangulation, these are DoFHandler and
MGDoFHandler. All of those define a set of iterators, allowing the user to
traverse the whole mesh, i.e. the set of cells, faces, edges, etc that
comprise the mesh, or portions of it. Therefore, they have somethings in
common. In fact, these iterators are instantiations or subclasses of the same
class TriaIterator (we do not include TriaRawIterator here, since it is only
for internal use).

The template signature of TriaIterator is
@code
  TriaIterator<dim, Accessor>
@endcode 

Usually, you will not use this definition directly, but employ one of the
typedefs below. Before going into this, let us first discuss the concept of
iterators, before delving into what the accessors do.

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
and decremented, but an operation like <tt>it=it+n</tt> takes a computing time
proportional to <tt>n</tt>, since it is implemented as a sequence of
<tt>n</tt> individual unit increments. Note that this is in contrast to the
next more specialized iterator concept, random access iterators, for which
access to an arbitrary object requires only constant time, rather than linear.


@section IteratorsAndSets Iterators as pointers into sets of objects

As mentioned above, iterators in deal.II can be considered as iterating over
all the objects that constitute a mesh. (These objects are lines, quads, and
hexes, and are each represented by a different kind of Accessor class; this
accessor is the second template argument in the code example above, and are
discussed in more detail below.) This suggests to view a triangulation as a
collection of cells and other objects that are held together by a certain data
structure that links all these objects, in the same was as a linked list is
the data structure that connects objects in a linear fashion.

Triangulations in deal.II can indeed be considered in this way. In particular,
they use the computational notion of a forest of regular trees to store their
data. This can be understood as follows: Consider the cells of the coarse mesh
as roots; then, if one of these coarse mesh cells is refined, it will have
2<sup>dim</sup> children, which in turn can, but do not have to have
2<sup>dim</sup> children of their own, and so on. This means, that each cell
of the coarse mesh can be considered the root of a binary tree (in qd), a
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
loop over all nodes of a forest, active iterators skip iterate over the
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
  KindIterator<dim,Accessor>
@endcode

Here, <tt>Kind</tt> determines what property an accessor needs to have to be
reached by this iterator (or omitted, for that matter). For example,
@code
  Iterator<dim,Accessor>
@endcode
iterates over all objects of kind Accessor that make up the mesh (for example
all cells, whether they are further refined and have children, or not), whereas
@code
  ActiveIterator<dim,Accessor>
@endcode
skips all objects that have children, i.e. objects that are not active.
Active iterators therefore operate on a strict subset of the objects
that normal iterators act on, namely those that possess the property that
they are active. Note that this is independent of the kind of object we
are operating on: all valid accessor classes have to provide the iterator
classes a method to find out whether they are active or not.

(For completeness, let us mention that there is a third kind of iterators: "raw
iterators" also traverse objects that are unused in the triangulation, but
allocated anyway for efficiency reasons. User code should never use raw
iterators, they are only for internal purposes of the library.)

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
  typename Trianguation<dim>::cell_iterator ti  = tria.active();
  typename DoFHandler<dim>::cell_iterator   di1 = dof1.active();
  typename DoFHandler<dim>::cell_iterator   di2 = dof2.active();
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
objects is the same for all iterators in the following snippet:
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
some useful object, in this case the Accessor. Accessing data that
characterizes a cell is always done through the Accessor, i.e. the expression
<tt>i-&gt;</tt> grants access to <b>all</b> attributes of this Accessor.

Examples of properties you can query from an iterator are
@code
  cell->vertex(1);
  line->child(0);
  hex->face(3);
  cell->at_boundary();
  face->boundary_indicator();
@endcode

Since dereferencing iterators yields accessor objects, these member functions
are declared in documented in the hierarchy of <code>TriaObjectAccessor</code>,
<code>CellAccessor</code>, <code>DoFObjectAccessor</code>,
<code>DoFCellAccessor</code>, <code>MGDoFObjectAccessor</code>, and
<code>MGDoFCellAccessor</code> classes. 


@section IteratorsTypedefs Iterators defined in the deal.II containers

Several classes in deal.II typedef iterators inside their class declarations.
The normal iterator types and calls to get them for cells and faces are:

<table border=1>
  <tr>
    <th>Container</th>
    <th>cell_iterator type</th>
    <th>function call</th>
  </tr>
  
  <tr>
    <th>Triangulation</th>
    <td>TriaIterator&lt;dim, CellAccessor&lt;dim&gt; &gt;</td>
    <td>triangulation.begin()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>TriaIterator&lt;dim, DoFCellAccessor&lt;dim&gt; &gt;</td>
    <td>dof_handler.begin()</td>
  </tr>

  <tr>
    <th>MGDoFHandler</th>
    <td>TriaIterator&lt;dim, MGDoFCellAccessor&lt;dim&gt; &gt;</td>
    <td>mg_dof_handler.begin()</td>
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
    <td>TriaIterator&lt;dim, TriaObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
    <td>triangulation.begin_face()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>TriaIterator&lt;dim, DoFObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
    <td>dof_handler.begin_face()</td>
  </tr>

  <tr>
    <th>MGDoFHandler</th>
    <td>TriaIterator&lt;dim, MGDoFObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
    <td>mg_dof_handler.begin_face()</td>
  </tr>
</table>


Likewise, active iterators are as follows:
<table border=1>
  <tr>
    <th>Container</th>
    <th>active_cell_iterator type</th>
    <th>function call</th>
  </tr>
  
  <tr>
    <th>Triangulation</th>
    <td>TriaActiveIterator&lt;dim, TriaCellAccessor&lt;dim&gt; &gt;</td>
    <td>triangulation.begin_active()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>TriaActiveIterator&lt;dim, DoFCellAccessor&lt;dim&gt; &gt;</td>
    <td>dof_handler.begin_active()</td>
  </tr>

  <tr>
    <th>MGDoFHandler</th>
    <td>TriaActiveIterator&lt;dim, MGDoFCellAccessor&lt;dim&gt; &gt;</td>
    <td>mg_dof_handler.begin_active()</td>
  </tr>
</table>


<table border=1>
  <tr>
    <th>Container</th>
    <th>active_face_iterator type</th>
    <th>function call</th>
  </tr>
  
  <tr>
    <th>Triangulation</th>
    <td>TriaActiveIterator&lt;dim, TriaObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
    <td>triangulation.begin_active_face()</td>
  </tr>

  <tr>
    <th>DoFHandler</th>
    <td>TriaActiveIterator&lt;dim, DoFObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
    <td>dof_handler.begin_active_face()</td>
  </tr>

  <tr>
    <th>MGDoFHandler</th>
    <td>TriaActiveIterator&lt;dim, MGDoFObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
    <td>mg_dof_handler.begin_active_face()</td>
  </tr>
</table>



*/
/**
@defgroup Accessors Accessor classes of the mesh iterators
*/

//@}

