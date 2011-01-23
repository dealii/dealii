//-----------------------  mapping_collection.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  mapping_collection.cc  ---------------------------

#include <base/memory_consumption.h>
#include <hp/mapping_collection.h>

DEAL_II_NAMESPACE_OPEN


namespace hp
{

  template<int dim, int spacedim>
  MappingCollection<dim,spacedim>::MappingCollection ()
  {}



  template<int dim, int spacedim>
  MappingCollection<dim,spacedim>::
  MappingCollection (const Mapping<dim,spacedim> &mapping)
  {
    mappings
      .push_back (std_cxx1x::shared_ptr<const Mapping<dim,spacedim> >(mapping.clone()));
  }



  template<int dim, int spacedim>
  MappingCollection<dim,spacedim>::
  MappingCollection (const MappingCollection<dim,spacedim> &mapping_collection)
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



  template<int dim, int spacedim>
  std::size_t
  MappingCollection<dim,spacedim>::memory_consumption () const
  {
    return (sizeof(*this) +
	    MemoryConsumption::memory_consumption (mappings));
  }



  template<int dim, int spacedim>
  void
  MappingCollection<dim,spacedim>::push_back (const Mapping<dim,spacedim> &new_mapping)
  {
    mappings
      .push_back (std_cxx1x::shared_ptr<const Mapping<dim,spacedim> >(new_mapping.clone()));
  }

//---------------------------------------------------------------------------



  template<int dim, int spacedim>
  MappingQ1<dim,spacedim>
  StaticMappingQ1<dim,spacedim>::mapping_q1;


  template<int dim, int spacedim>
  MappingCollection<dim,spacedim>
  StaticMappingQ1<dim,spacedim>::mapping_collection
  = MappingCollection<dim,spacedim>(mapping_q1);

}



// explicit instantiations
#include "mapping_collection.inst"


DEAL_II_NAMESPACE_CLOSE
