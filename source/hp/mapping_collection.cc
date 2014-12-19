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
#include <deal.II/hp/mapping_collection.h>

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
    .push_back (std_cxx11::shared_ptr<const Mapping<dim,spacedim> >(mapping.clone()));
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
    .push_back (std_cxx11::shared_ptr<const Mapping<dim,spacedim> >(new_mapping.clone()));
  }

//---------------------------------------------------------------------------


  namespace
  {
    /**
     * Create and return a reference to a static MappingQ1 object. We can't
     * use the one in ::StaticMappingQ1 to initialize the static object below
     * since we can't make sure that the constructor for that object is run
     * before we want to use the object (when constructing mapping_collection
     * below).  Therefore we create a helper function which returns a
     * reference to a static object that will be constructed the first time
     * this function is called.
     */
    template<int dim, int spacedim>
    MappingQ1<dim,spacedim> &
    get_static_mapping_q1()
    {
      static MappingQ1<dim,spacedim> mapping;
      return mapping;
    }
  }

  template<int dim, int spacedim>
  MappingCollection<dim,spacedim>
  StaticMappingQ1<dim,spacedim>::mapping_collection
    = MappingCollection<dim,spacedim>(get_static_mapping_q1<dim,spacedim>());

}



// explicit instantiations
#include "mapping_collection.inst"


DEAL_II_NAMESPACE_CLOSE
