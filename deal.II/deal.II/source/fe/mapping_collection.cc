//-----------------------  mapping_collection.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  mapping_collection.cc  ---------------------------

#include <base/memory_consumption.h>
#include <fe/mapping_collection.h>


namespace hp
{

  template <int dim>
  MappingCollection<dim>::MappingCollection () :
                  single_mapping (false)
  {
  }


  template <int dim>
  MappingCollection<dim>::MappingCollection (const Mapping<dim> &mapping) :
                  single_mapping (true)
  {
    mappings
      .push_back (boost::shared_ptr<const Mapping<dim> >(mapping.clone()));
  }


  template <int dim>
  inline
  const Mapping<dim> &
  MappingCollection<dim>::operator[] (const unsigned int index) const
  {
    if (single_mapping)
      return *mappings[0];
    else
      {
	Assert (index < mappings.size (),
		ExcIndexRange (index, 0, mappings.size ()));
	return *mappings[index];
      }
  }


  template <int dim>
  unsigned int
  MappingCollection<dim>::memory_consumption () const
  {
    return (sizeof(*this) +
	    MemoryConsumption::memory_consumption (mappings));
  }


  template <int dim>
  unsigned int MappingCollection<dim>::
  add_mapping (const Mapping<dim> &new_mapping)
  {
                                     // A MappingCollection, which was
                                     // initialized as single
                                     // MappingCollection cannot administrate
                                     // other Mappings.
    Assert (!single_mapping, ExcNotInitialized ());
    
    mappings
      .push_back (boost::shared_ptr<const Mapping<dim> >(new_mapping.clone()));

    return mappings.size ();
  }


// explicit instantiations
  template class MappingCollection<deal_II_dimension>;

  
}
