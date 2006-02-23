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
#include <fe/mapping_q1.h>
#include <fe/fe.h>

#include <vector>
#include <boost/shared_ptr.hpp>


//TODO[WB]: The conversion constructor should really be 'explicit'

namespace hp
{
/**
 * This class implements a collection of mapping objects in the same way as
 * the hp::FECollection implements a collection of finite element classes.
 *
 * It implements the concepts stated in the @ref hpcollection module described
 * in the doxygen documentation.
 * 
 * Although it is recommended to supply an appropriate mapping for
 * each finite element kind used in a hp-computation, the
 * MappingCollection class implements a conversion constructor from a
 * single mapping.  Therefore it is possible to offer only a single
 * mapping to the hp::FEValues class instead of a
 * hp::MappingCollection. This is for the convenience of the user, as
 * many simple geometries do not require different mappings along the
 * boundary to achieve optimal convergence rates.  Hence providing a
 * single mapping object will usually suffice. See the hp::FEValues
 * class for the rules which mapping will be selected for a given
 * cell.
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
                                        * Default constructor. Leads
                                        * to an empty collection that
                                        * can later be filled using
                                        * push_back().
                                        */
      MappingCollection ();

                                       /**
                                        * Conversion constructor. This
                                        * constructor creates a
                                        * MappingCollection from a
                                        * single mapping. More
                                        * mappings can be added with
                                        * push_back(), if desired,
                                        * though it would probably be
                                        * clearer to add all mappings
                                        * the same way.
                                        */
      MappingCollection (const Mapping<dim> &mapping);

                                       /**
                                        * Copy constructor.
                                        */
      MappingCollection (const MappingCollection<dim> &mapping_collection);
      
                                       /**
                                        * Adds a new mapping to the
                                        * MappingCollection.  The
                                        * mappings have to be added in
                                        * the order of the
                                        * active_fe_indices. Thus the
                                        * reference to the mapping
                                        * object for active_fe_index 0
                                        * has to be added first,
                                        * followed by the mapping
                                        * object for active_fe_index
                                        * 1.
                                        */
      void push_back (const Mapping<dim> &new_mapping);

                                       /**
                                        * Returns the mapping object
                                        * which was specified by the
                                        * user for the active_fe_index
                                        * which is provided as a
                                        * parameter to this method.
                                        *
                                        * @pre @p index must be
                                        * between zero and the number
                                        * of elements of the
                                        * collection.
                                        */
      const Mapping<dim> &
      operator[] (const unsigned int index) const;

                                       /**
                                        * Returns the number of
                                        * mapping objects stored in
                                        * this container.
                                        */
      unsigned int size () const;
    
                                       /**
                                        * Determine an estimate for the
                                        * memory consumption (in bytes)
                                        * of this object.
                                        */
      unsigned int memory_consumption () const;

    private:
                                       /**
                                        * The real container, which stores
                                        * pointers to the different Mapping
                                        * objects.
                                        */
      std::vector<boost::shared_ptr<const Mapping<dim> > > mappings;
  };


/**
 * In order to avoid creation of static MappingQ1 objects at several
 * places in the library (in particular in backward compatibility
 * functions), we define a static MappingQ1 objects once and for all
 * places where it is needed.
 */
  template <int dim>
  struct StaticMappingQ1
  {
      static MappingCollection<dim> mapping_collection;
  };


/* --------------- inline functions ------------------- */

  template <int dim>
  inline
  unsigned int
  MappingCollection<dim>::size () const 
  {
    return mappings.size();
  }



  template <int dim>
  inline
  const Mapping<dim> &
  MappingCollection<dim>::operator[] (const unsigned int index) const
  {
    Assert (index < mappings.size (),
            ExcIndexRange (index, 0, mappings.size ()));
    return *mappings[index];
  }
  
} // namespace hp


#endif
