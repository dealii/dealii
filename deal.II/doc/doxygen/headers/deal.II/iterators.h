/**
@page Iterators Iterators on mesh like containers

deal.II has several classes which are understood conceptionally as
meshes. Apart from the obvious Triangulation, these are DoFHandler and
MGDoFHandler. All of those define a set of iterators, allowing the
user to traverse the whole mesh or portions of it. Therefore, they
have somethings in common. In fact, these iterators are instantiations
or subclasses of the same class TriaIterator (we do not include
TriaRawIterator here, since it is only for internal use),

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
the DoFHandlers are handling different finite elements. This
distinction is only in the Accessor.

The standard loops are
<dl>
<dt>TriaIterator</tt>
<dd>Traverse all cells on all levels</dd>

<dt>TriaActiveIterator</tt>
<dd>Loop over @ref GlossActive "active cells" only</dd> 
</dl>

@subsection IteratorsAccessors Accessors

We just saw that iterators are not very clever and just know how to
march through a mesh. Whatever you actually can do with this iterator
is determined by the Accessor. The magic is, that for any iterator
<tt>i</tt>, the term <tt>i-&gt;</tt> grants access to <b>all</b>
attributes of this Accessor.


@section IteratorsTypedefs Iterators defined in the containers

The standard iterators are typedefed inside the classes. These are

<table border=1>
<tr><th></th>
<th>cell_iterator</th>
<th>face_iterator</th>
</table>
