//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**

@page UpdateFlagsEssay The interplay of UpdateFlags, Mapping and FiniteElement in FEValues

@section update_flags_Intro Introduction

The class FEValues was introduced because most of the time values of
shape functions etc. are used in Quadrature points only and their
computation might be sped up by doing it all at once. Furthermore,
some quantities, like the shape function values of standard scalar
finite elements do not depend on the actual mesh cell at all and do
not have to be reevaluated on each cell. The interface between
FEValues and FiniteElement should be aware of such things and do the
best thing possible.

Now, in order not to compute too much, FEValues has to be told what to
update. This is where UpdateFlags comes into play. And each FEValues
uses two sets of UpdateFlags internally: one that describes the values
that can be computed on the reference cell only and one requiring
recomputation on each cell. Furthermore, since the actual computations
are performed by Mapping and FiniteElement objects, both of them have
UpdateFlags.

@section update_flags_Once Update once or each

Your desired set of values for computation given, the FiniteElement
and Mapping involved have to perform the following tasks:
<ol>
<li> Are any additional values required in order to compute these
values? If so, add these flags to the current set. For instance, the
derivative od a standard scalar element requires the inverse of the
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

@section update_flags_Functions Generation of the actual data

The computation of the fields accessible through FEValues proceeds in
two different steps:

@subsection update_flags_Init Initialization

Functions of FEValues classes involved:
<ul>
<li>FEValues::initialize()
<li>FEFaceValues::initialize()
<li>FESubfaceValues::initialize()
</ul>

These function is called by the constructor of FEValues, FEFaceValues
and FESubfaceValues, respectively, and allows the
Mapping and FiniteElement objects to set up internal data
structures. These structures are internal in the sense that FEValues
only forwards them to Mapping and FiniteElement if necessary, but does
not access their contents. FEValues just stores a pointer to their
base class Mapping::InternalDataBase. Concrete Mapping and
FiniteElement classes will derive their own sutructures containing all
the data needed to speed up updating the fields of FEValues on a cell
or face.

The functions actually providing these internal data structures are:
<ul>
<li>Mapping::get_data()
<li>Mapping::get_face_data()
<li>Mapping::get_subface_data()
<li>FiniteElement::get_data()
<li>FiniteElement::get_face_data()
<li>FiniteElement::get_subface_data()
</ul>

@subsection update_flags_Reinit Reinitialization for a mesh cell

Functions of FEValues classes involved:
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
