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
  MappingCollection<dim>::MappingCollection ()
  {}


  
  template <int dim>
  MappingCollection<dim>::
  MappingCollection (const Mapping<dim> &mapping)
  {
    mappings
      .push_back (boost::shared_ptr<const Mapping<dim> >(mapping.clone()));
  }



  template <int dim>
  MappingCollection<dim>::
  MappingCollection (const MappingCollection<dim> &mapping_collection)
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
                  mappings (mapping_collection.mappings)
  {}
  
  

  template <int dim>
  unsigned int
  MappingCollection<dim>::memory_consumption () const
  {
    return (sizeof(*this) +
	    MemoryConsumption::memory_consumption (mappings));
  }

  

  template <int dim>
  void
  MappingCollection<dim>::push_back (const Mapping<dim> &new_mapping)
  {
    mappings
      .push_back (boost::shared_ptr<const Mapping<dim> >(new_mapping.clone()));
  }

//---------------------------------------------------------------------------


  template <int dim>
  MappingCollection<dim>
  StaticMappingQ1<dim>::mapping_collection
  = MappingCollection<dim>(::StaticMappingQ1<dim>::mapping);

// explicit instantiations
  template class MappingCollection<deal_II_dimension>;
  template struct StaticMappingQ1<deal_II_dimension>;
  
}
