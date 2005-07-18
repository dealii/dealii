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

Usually, you will not use this definition directly, but employ one of
the typedefs below.

@section IteratorsDifferences Distinguishing between iterators

The iterators discussed are all of the form
@code
  LoopIterator<dim, Accessor>
@endcode

Here, <tt>Loop</tt> determines, which cells are reached (or omitted,
for that matter). This means, independent of the accessor type, this
part of the definition of the iterator determines the meaning of the
increment operator.

@subsection IteratorsLoops The action of the iterator itself

All iterators with the same <tt>Loop</tt> and iterating over the
same kind of geometrical objects traverse the mesh in the same
order. Take this code example:
@code
  Triangulation<dim> tria;
  DoFHandler<dim>    dof1(tria);
  DoFHandler<dim>    dof2(tria);
  ...
  Trianguation<dim>::active_cell_iterator ti  = tria.begin_active();
  DoFHandler<dim>::active_cell_iterator   di1 = dof1.begin_active();
  DoFHandler<dim>::active_cell_iterator   di2 = dof2.begin_active();
  ...
  while (ti != tria.end())
  {
    do_something_with_iterators(ti, di1, di2);
    ++ti;
    ++di1;
    ++di2;
  }
@endcode

Here, all iterators will always point to the same mesh cell, even if
the DoFHandlers are handling different finite elements: they all access cells in the same order, the difference is only in the Accessor.

The standard loops are
<dl>
<dt>TriaIterator</dt>
<dd>Traverse all cells on all levels</dd>

<dt>TriaActiveIterator</dt>
<dd>Loop over @ref GlossActive "active cells" only</dd> 
</dl>

In addition, there are "raw iterators" that also traverse objects that are
unused in the triangulation, but allocated anyway for efficiency reasons. User
code should never use raw iterators, they are only for internal purposes of
the library.


@subsection IteratorsAccessors Accessors

Iterators are like pointers: they can be incremented and decremented, but they
are really rather dumb. Their magic only lies in the fact that they point to
some useful object, in this case the Accessor. Accessing data that
characterizes a cell is always done through the Accessor, i.e. the expression
<tt>i-&gt;</tt> grants access to <b>all</b> attributes of this Accessor.


@section IteratorsTypedefs Iterators defined in the containers

The standard iterators are typedefed inside the classes. These are

<table border=1>
<tr><th></th>
<th>cell_iterator</th>
<th>face_iterator</th>
</tr>
<tr>
<th>Triangulation</th>
<td>TriaIterator&lt;dim, CellAccessor&lt;dim&gt; &gt;</td>
<td>TriaIterator&lt;dim, TriaObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
</tr>
<tr>
<th>DoFHandler</th>
<td>TriaIterator&lt;dim, DoFCellAccessor&lt;dim&gt; &gt;</td>
<td>TriaIterator&lt;dim, DoFObjectAccessor&lt;dim-1, dim&gt; &gt;</td>
</tr>
</table>

@section IteratorsFiltered Filtered iterators

Filtered iterators restrict the scope of existing iterators even
more. For instance, you could imagine to iterate over the the subset
of those @ref GlossActive "active cells" having their user flag set or
belonging to a certain subdomain.

This is achieved by using an object of type FilteredIterator
&lt;BaseIterator&gt;, where BaseIterator usually is one of the
standard iterators discussed above.

The FilteredIterator gets an additional Predicate in its constructor
and will skip all objects where this Predicate evaluates to
<tt>false</tt>. A collection of already implemented predicates can be
found in the namespace IteratorFilters.


@section IteratorsAndSets Iterators as pointers into sets of objects

As mentioned above, iterators in deal.II can be considered as iterating over
all the cells (or faces, lines, etc) that constitute a mesh. This suggests to
view a triangulation as a collection of cells and other objects that are held
together by a certain data structure that links all these objects, in the same
was as a linked list is the data structure that connects objects in a linear
fashion.

Triangulations in deal.II can indeed be considered in this way. In particular,
they use the computational notion of a forrest of regular trees to store their
data. This can be understood as follows: Consider the cells of the coarse mesh
as roots; then, if one of these coarse mesh cells is refined, it will have
2<sup>dim</sup> children, which in turn can, but do not have to have
2<sup>dim</sup> children of their own, and so on. This means, that each cell
of the coarse mesh can be considered the root of a binary tree (in qd), a
quadtree (in 2d), or an octree (in 3d). The collection of these trees
emanating from the cells of the coarse mesh then constitutes the forrest that
completely describes the triangulation, including all of its active and
inactive cells. In particular, the active cells are those terminal nodes in
the tree that have no decendants, i.e. cells which are not further
refined. Correspondingly, inactive cells correspond to nodes in the tree with
descendents, i.e. cells that are further refined.

A triangulation contains forrests for lines (each of which may have 2
children), quads (each with possibly four children), and hexes (each with no
or 8 children). Depending on the dimension, these objects are also termed
cells or faces.

Iterators loop over the elements of such forrests. While the usual iterators
loop over all nodes of a forrest, active iterators skip iterate over the
elements in the same order, but skip all non-active entries and therefore only
visit terminal nodes (i.e. active cells, faces, etc). There are many ways to
traverse the elements of a forrest, for example breadth first or depth
first. Depending on the type of data structure used to store the forrest, some
ways are more efficient than others. At present, the way iterators traverse
forrests in deal.II is breadth first. I.e., iterators first visit all the
elements (cells, faces, etc) of the coarse mesh before moving on to all the
elements of the immediate level, i.e. the immediate children of the coarse
mesh objects; after this come the grandchildren of the coarse mesh, and so on.
However, it must be noted that programs should not rely on this particular
order of traversing a tree: this is considered an implementation detail that
can change between versions, even if we consider this an unlikely option at
the present time.

*/
/**
@defgroup Accessors Accessor classes of the mesh iterators
*/

//@}

