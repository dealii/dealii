// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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

@page UpdateFlagsEssay The interplay of UpdateFlags, Mapping and FiniteElement in FEValues

<h2>Introduction</h2>

In order to compute local contributions of an individual to the global
matrix and right hand side, integrals are usually transformed to the
reference cell. For example, for the Laplace matrix, on each cell we
have to compute
@f[
  A^K_{ij} = \sum_{q}J^{-1}(\hat{\bf x}_q) \hat \nabla \varphi_i(\hat{\bf x}_q) \cdot
  J^{-1}(\hat{\bf x}_q) \hat \nabla \varphi_j(\hat{\bf x}_q)\ |\textrm{det}\ J(\hat{\bf x}_q)|
  w_q,
@f]
where a hat indicates reference coordinates, and $J(\hat{\bf
x}_q)$ is the Jacobian
$\frac{\partial F_K(\hat{\bf x}_q)}{\partial\bf \hat x}$ of the mapping,
evaluated at a quadrature point $\hat{\bf x}_q$ on the reference cell.

In order to evaluate such an expression in an application code, we
have to access three different kinds of objects: a quadrature
object that describes locations $\hat{\bf x}_q$ and weights $w_q$ of
quadrature points on the reference cell; a finite element object that
describes the gradients $\hat\nabla \varphi_i(\hat{\bf x}_q)$ of shape
functions on the unit cell; and a mapping object that provides the
Jacobian as well as its determinant. Dealing with all these
objects would be cumbersome and error prone.

On the other hand, these three kinds of objects almost always appear together,
and it is in fact very rare for deal.II application codes to do anything with
quadrature, finite element, or mapping objects besides using them together.
For this reason, we have introduced the FEValues abstraction
combining information on the shape functions, the geometry of the actual mesh
cell and a quadrature rule on a reference cell. Upon construction it takes one
object of each of the three mentioned categories. Later, it can be
"re-initialized" for a concrete grid cell and then provides mapped quadrature
points and weights, mapped shape function values and derivatives as well as
some properties of the transformation from the reference cell to the actual
mesh cell.

Since computation of any of these values is potentially expensive (for
example when using high order mappings with high order elements), the
FEValues class only computes what is explicitly asked for. To this
end, it takes a list of flags of type UpdateFlags at construction time
specifying which quantities should be updated each time a cell is
visited. In addition, allowing further optimizations, the functions
filling the data fields of FEValues are able to distinguish between
values that have to be recomputed on each cell (for example mapped
gradients) and quantities that do not change from cell to cell (for
example the values of shape functions of the usual $Q_k$
finite elements at the same quadrature points on different cells; this
property does not hold for the shape functions of Raviart-Thomas
elements, however, which must be rotated with the local cell).

At construction time, each FEValues object splits the given
UpdateFlags into two sets: those flags that describe quantities that
can be computed on the reference cell once at the beginning
("update_once"), and those that require recomputation on each cell
("update_each").

Furthermore, each FEValues object is associated with a Mapping and a
FiniteElement object that do the actual computations. It therefore
keeps a set up UpdateFlags for both update_once and update_each
operations for both the Mapping and the FiniteElement object it is
associated with.


<h2>Update once or each</h2>

Sometimes in order to compute one quantity, something else also has to be
computed. For example, if only update_values is requested by the user of
a FEValues object and if the associated finite element is of type FE_Q,
then <code>update_once=update_values</code> and <code>update_each=0</code>
for the FiniteElement computations and <code>update_once=update_each=0</code>
for the Mapping object: we can compute the values of the shape function
at the quadrature points of each cell once on the reference cell, and they
will have the same value at the quadrature point of a real cell. There is
nothing the finite element object has to do when we move to the next cell
(these would be update_each calculations) and there is nothing the mapping
has to do, either up front (update_once) or on each cell (update_each).

However, this is not the case if we used a FE_RaviartThomas element: there,
computing the values of the shape functions on a cell involves knowing the
Jacobian of the mapping, and so if a user requests <code>update_values</code>
of a FEValues object that uses a FE_RaviartThomas element, then we can set
<code>update_once=update_values</code> and <code>update_each=0</code>
for the FiniteElement, but need to set <code>update_once=0</code>
<code>update_each=update_jacobians</code> for the Mapping object.

To accommodate this structure, at the time a FEValues object is constructed,
it asks both the FiniteElement and the Mapping object it uses the following:
<ol>
<li> Are any additional values required in order to compute the
currently required values? If so, add these flags to the current set.
Another example to the one above would be that the
derivative of a standard scalar element requires the inverse of the
Jacobian of the Mapping.
<li> Given the enhanced set, determine the subsets of values that is
performed on the reference cell only and on each cell, respectively.
</ol>

In order to compute this, there are functions Mapping::update_once()
and Mapping::update_each() as well as FiniteElement::update_once() and
FiniteElement::update_each(). All of them accept UpdateFlags as the
set of desired flags and return UpdateFlags as the set of required
flags on the reference and on each cell, respectively. Additionally,
FiniteElement::update_once() should set all flags for values that are
required from the Mapping. The function in FEValues computing the
actual set of flags from the desired one looks like this:
<code>
  flags |= fe->update_once (flags)
        |  fe->update_each (flags);
  flags |= mapping->update_once (flags)
        |  mapping->update_each (flags);
</code>
That is, a FiniteElement can set additional flags which are honored by
the Mapping.



<h2>Generation of the actual data</h2>

As outlined above, data is computed at two different times: once at
the beginning on the reference cell, and once whenever we move to an
actual cell. The functions involved in each of these steps are
discussed next:


<h3>Initialization</h3>

Computing data on the reference cell before we even visit the first
real cell is a two-step process. First, the constructor of FEValues,
FEFaceValues and FESubfaceValues, respectively, need to allow the
Mapping and FiniteElement objects to set up internal data
structures. These structures are internal in the following sense: the
FEValues object asks the finite element and mapping objects to create
an object of type FiniteElement::InternalDataBase and
Mapping::InternalDataBase each; the actual finite element and mapping
class may in fact create objects of a derived type if they wish to
store some data beyond what these base classes already provide. The
functions involved in this are
<ul>
<li>Mapping::get_data()
<li>Mapping::get_face_data()
<li>Mapping::get_subface_data()
<li>FiniteElement::get_data()
<li>FiniteElement::get_face_data()
<li>FiniteElement::get_subface_data()
</ul>

The FEValues object then takes over ownership of these objects and will
destroy them at the end of the FEValues object's lifetime. After this,
the FEValues object asks the FiniteElement and Mapping objects to fill
these InternalDataBase objects with the data that pertains to what
can and needs to be computed on the reference cell. This is done in these
functions:
<ul>
<li>FEValues::initialize()
<li>FEFaceValues::initialize()
<li>FESubfaceValues::initialize()
</ul>


<h3>Reinitialization for a mesh cell</h3>

Once initialization is over and we call FEValues::reinit, FEFaceValues::reinit
or FESubfaceValues::reinit to move to a concrete cell or face, we need
to calculate the update_each kinds of data. This done in the following
functions:
<ul>
<li>FEValues::reinit() calls Mapping::fill_fe_values(), then FiniteElement::fill_fe_values()
<li>FEFaceValues::reinit() calls Mapping::fill_fe_face_values(), then FiniteElement::fill_fe_face_values()
<li>FESubfaceValues::reinit() calls Mapping::fill_fe_subface_values(),
thenFiniteElement::fill_fe_subface_values()
</ul>

This is, where the actual data fields for FEValues, stored in
FEValuesData objects is computed. These functions call the function in
Mapping first, such that all the mapping data required by the finite
element is available. Then, the FiniteElement function is called.

When this happens for the first time after initialization, all the
values specified by Mapping::InternalDataBase::update_once or
Mapping::InternalDataBase::update_each are filled. After that, only
the values specified by Mapping::InternalDataBase::update_each will be
updated.

*/
