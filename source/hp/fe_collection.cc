// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


#include <deal.II/base/memory_consumption.h>
#include <deal.II/hp/fe_collection.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  template <int dim, int spacedim>
  FECollection<dim,spacedim>::FECollection ()
  {}



  template <int dim, int spacedim>
  FECollection<dim,spacedim>::FECollection (const FiniteElement<dim,spacedim> &fe)
  {
    push_back (fe);
  }



  template <int dim, int spacedim>
  FECollection<dim,spacedim>::
  FECollection (const FECollection<dim,spacedim> &fe_collection)
    :
    Subscriptor (),
    // copy the array
    // of shared
    // pointers. nothing
    // bad should
    // happen -- they
    // simply all point
    // to the same
    // objects, and the
    // last one to die
    // will delete the
    // mappings
    finite_elements (fe_collection.finite_elements)
  {}



  template <int dim, int spacedim>
  void FECollection<dim,spacedim>::push_back (const FiniteElement<dim,spacedim> &new_fe)
  {
    // check that the new element has the right
    // number of components. only check with
    // the first element, since all the other
    // elements have already passed the test
    // against the first element
    if (finite_elements.size() != 0)
      Assert (new_fe.n_components() == finite_elements[0]->n_components(),
              ExcMessage ("All elements inside a collection need to have the "
                          "same number of vector components!"));

    finite_elements
    .push_back (std_cxx11::shared_ptr<const FiniteElement<dim,spacedim> >(new_fe.clone()));
  }



  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim,spacedim>::
  component_mask (const FEValuesExtractors::Scalar &scalar) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].component_mask(scalar),
              ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim,spacedim>::
  component_mask (const FEValuesExtractors::Vector &vector) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].component_mask(vector),
              ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim,spacedim>::
  component_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].component_mask(sym_tensor),
              ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim,spacedim>::
  component_mask (const BlockMask &block_mask) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(block_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].component_mask(block_mask),
              ExcMessage ("Not all elements of this collection agree on what "
                          "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim,spacedim>::
  block_mask (const FEValuesExtractors::Scalar &scalar) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].block_mask(scalar),
              ExcMessage ("Not all elements of this collection agree on what "
                          "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim,spacedim>::
  block_mask (const FEValuesExtractors::Vector &vector) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].block_mask(vector),
              ExcMessage ("Not all elements of this collection agree on what "
                          "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim,spacedim>::
  block_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].block_mask(sym_tensor),
              ExcMessage ("Not all elements of this collection agree on what "
                          "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  BlockMask
  FECollection<dim,spacedim>::
  block_mask (const ComponentMask &component_mask) const
  {
    Assert (size() > 0,
            ExcMessage ("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(component_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c=1; c<size(); ++c)
      Assert (mask == (*this)[c].block_mask(component_mask),
              ExcMessage ("Not all elements of this collection agree on what "
                          "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::n_blocks () const
  {
    Assert (finite_elements.size () > 0, ExcNoFiniteElements());

    const unsigned int nb = finite_elements[0]->n_blocks ();
    for (unsigned int i=1; i<finite_elements.size(); ++i)
      Assert (finite_elements[i]->n_blocks() == nb,
              ExcMessage ("Not all finite elements in this collection have "
                          "the same number of components."));

    return nb;
  }



  template <int dim, int spacedim>
  std::size_t
  FECollection<dim,spacedim>::memory_consumption () const
  {
    std::size_t mem
      = (sizeof(*this) +
         MemoryConsumption::memory_consumption (finite_elements));
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      mem += finite_elements[i]->memory_consumption();

    return mem;
  }
}



// explicit instantiations
#include "fe_collection.inst"


DEAL_II_NAMESPACE_CLOSE
