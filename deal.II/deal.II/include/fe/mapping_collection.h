//------------------------  mapping_collection.h  --------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------  mapping_collection.h  --------------------------
#ifndef __deal2__mapping_collection_h
#define __deal2__mapping_collection_h

#include <base/config.h>
#include <base/subscriptor.h>
#include <fe/mapping.h>
#include <fe/fe.h>

#include <vector>
#include <boost/shared_ptr.hpp>


namespace hp
{
/**
 * This class implements a collection of mapping objects in the same way as
 * the hp::FECollection implements a collection of finite element classes.
 *
 * It implements the concepts stated in the @ref hpcollection module described
 * in the doxygen documentation.
 * 
 * Although it is recommended to supply an appropriate mapping for each finite
 * element kind used in a hp-computation, the MappingCollection class
 * implements a conversion constructor from a single mapping.  Therefore it is
 * possible to offer only a single mapping to the hp::FEValues class instead of
 * a hp::MappingCollection. This is for the convenience of the user, as many
 * simple geometries do not require different mappings along the boundary to
 * achieve optimal convergence rates.  Hence providing a single mapping object
 * will usually suffice.
 * 
 * @ingroup hp hpcollection
 * 
 * @author Oliver Kayser-Herold, 2005
 */
  template <int dim>
  class MappingCollection : public Subscriptor
  {
    public:
                                       /**
                                        * Default constructor. Initialises
                                        * this MappingCollection.
                                        */
      MappingCollection ();

                                       /**
                                        * Conversion constructor. This
                                        * constructor creates a MappingCollection
                                        * from a single mapping. In
                                        * the newly created MappingCollection, this
                                        * mapping is used for all active_fe
                                        * indices.
                                        */
      MappingCollection (const Mapping<dim> &mapping);

                                       /**
                                        * Returns the number of mapping
                                        * objects stored in this container.
                                        */
      unsigned int size () const;

                                       /**
                                        * Returns the mapping object which
                                        * was specified by the user for the
                                        * active_fe_index which is provided
                                        * as a parameter to this method.
                                        *
                                        * @pre @p index must be between zero
                                        * and the number of elements of the
                                        * collection.
                                        */
      const Mapping<dim> &
      operator[] (const unsigned int index) const;
    
                                       /**
                                        * Determine an estimate for the
                                        * memory consumption (in bytes)
                                        * of this object.
                                        */
      unsigned int memory_consumption () const;

                                       /**
                                        * Adds a new mapping to the
                                        * MappingCollection.  The mappings
                                        * have to be added in the order of the
                                        * active_fe_indices. Thus the
                                        * reference to the mapping object for
                                        * active_fe_index 0 has to be added
                                        * first, followed by the mapping
                                        * object for active_fe_index 1.
                                        */
      void push_back (const Mapping<dim> &new_mapping);

    private:
                                       /**
                                        * Upon construction of a
                                        * <tt>MappingCollection</tt> the later
                                        * functionality of the class is
                                        * specified.  Either it is a real
                                        * collection, which provides different
                                        * mappings for each active_fe_index or
                                        * its a "unreal" collection, which
                                        * returns a the same mapping for all
                                        * active_fe_indices.  This boolean
                                        * remembers which type this object is.
                                        */
      const bool single_mapping;

                                       /**
                                        * The real container, which stores
                                        * pointers to the different Mapping
                                        * objects.
                                        */
      std::vector<boost::shared_ptr<const Mapping<dim> > > mappings;
  };



/* --------------- inline functions ------------------- */

  template <int dim>
  inline
  unsigned int
  MappingCollection<dim>::size () const 
  {
    return mappings.size();
  }
  
} // namespace hp


#endif
