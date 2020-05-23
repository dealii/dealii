// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_mapping_collection_h
#define dealii_mapping_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * This class implements a collection of mapping objects in the same way as
   * the hp::FECollection implements a collection of finite element classes.
   *
   * It implements the concepts stated in the
   * @ref hpcollection
   * module described in the doxygen documentation.
   *
   * Although it is recommended to supply an appropriate mapping for each
   * finite element kind used in a hp-computation, the MappingCollection class
   * implements a conversion constructor from a single mapping.  Therefore it
   * is possible to offer only a single mapping to the hp::FEValues class
   * instead of a hp::MappingCollection. This is for the convenience of the
   * user, as many simple geometries do not require different mappings along
   * the boundary to achieve optimal convergence rates.  Hence providing a
   * single mapping object will usually suffice. See the hp::FEValues class
   * for the rules which mapping will be selected for a given cell.
   *
   * @ingroup hp hpcollection
   *
   * @author Oliver Kayser-Herold, 2005
   */
  template <int dim, int spacedim = dim>
  class MappingCollection : public Subscriptor
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    MappingCollection() = default;

    /**
     * Conversion constructor. This constructor creates a MappingCollection
     * from a single mapping. More mappings can be added with push_back(), if
     * desired, though it would probably be clearer to add all mappings the
     * same way.
     */
    explicit MappingCollection(const Mapping<dim, spacedim> &mapping);

    /**
     * Copy constructor.
     */
    MappingCollection(
      const MappingCollection<dim, spacedim> &mapping_collection);

    /**
     * Add a new mapping to the MappingCollection. Generally, you will
     * want to use the same order for mappings as for the elements of
     * the hp::FECollection object you use. However, the same
     * considerations as discussed with the hp::QCollection::push_back()
     * function also apply in the current context.
     *
     * This class creates a copy of the given mapping object, i.e., you can
     * do things like <tt>push_back(MappingQ<dim>(3));</tt>. The internal copy
     * is later destroyed by this object upon destruction of the entire
     * collection.
     */
    void
    push_back(const Mapping<dim, spacedim> &new_mapping);

    /**
     * Return the mapping object which was specified by the user for the
     * active_fe_index which is provided as a parameter to this method.
     *
     * @pre @p index must be between zero and the number of elements of the
     * collection.
     */
    const Mapping<dim, spacedim> &operator[](const unsigned int index) const;

    /**
     * Return the number of mapping objects stored in this container.
     */
    unsigned int
    size() const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

  private:
    /**
     * The real container, which stores pointers to the different Mapping
     * objects.
     */
    std::vector<std::shared_ptr<const Mapping<dim, spacedim>>> mappings;
  };


  /**
   * Many places in the library by default use (bi-,tri-)linear mappings
   * unless users explicitly provide a different mapping to use. In these
   * cases, the called function has to create a $Q_1$ mapping object, i.e., an
   * object of kind MappingQGeneric(1). This is costly. It would also be
   * costly to create such objects as static objects in the affected
   * functions, because static objects are never destroyed throughout the
   * lifetime of a program, even though they only have to be created once the
   * first time code runs through a particular function.
   *
   * In order to avoid creation of (static or dynamic) $Q_1$ mapping objects
   * in these contexts throughout the library, this class defines a static
   * collection of mappings with a single $Q_1$ mapping object. This
   * collection can then be used in all of those places where such a
   * collection is needed.
   */
  template <int dim, int spacedim = dim>
  struct StaticMappingQ1
  {
  public:
    /**
     * The publicly available static $Q_1$ mapping collection object.
     */
    static MappingCollection<dim, spacedim> mapping_collection;
  };


  /* --------------- inline functions ------------------- */

  template <int dim, int spacedim>
  inline unsigned int
  MappingCollection<dim, spacedim>::size() const
  {
    return mappings.size();
  }



  template <int dim, int spacedim>
  inline const Mapping<dim, spacedim> &MappingCollection<dim, spacedim>::
                                       operator[](const unsigned int index) const
  {
    AssertIndexRange(index, mappings.size());
    return *mappings[index];
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
