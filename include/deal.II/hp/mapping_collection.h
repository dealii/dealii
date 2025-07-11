// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_collection_h
#define dealii_mapping_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/collection.h>

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
   * topic described in the doxygen documentation.
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
   */
  template <int dim, int spacedim = dim>
  class MappingCollection : public Collection<Mapping<dim, spacedim>>
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
     * Constructor. This constructor creates a MappingCollection from one or
     * more mapping objects passed to the constructor. For this
     * call to be valid, all arguments need to be of types derived
     * from class Mapping<dim,spacedim>.
     */
    template <class... MappingTypes>
    explicit MappingCollection(const MappingTypes &...mappings);

    /**
     * Copy constructor.
     */
    MappingCollection(
      const MappingCollection<dim, spacedim> &mapping_collection);

    /**
     * Move constructor.
     *
     * @note The implementation of standard datatypes may change with different
     * libraries, so their move members may or may not be flagged non-throwing.
     * We need to explicitly set the noexcept specifier according to its
     * member variables to still get the performance benefits (and to satisfy
     * clang-tidy).
     */
    MappingCollection(MappingCollection<dim, spacedim> &&) noexcept(
      std::is_nothrow_move_constructible_v<
        std::vector<std::shared_ptr<const Mapping<dim, spacedim>>>>
        &&std::is_nothrow_move_constructible_v<std::function<
          unsigned int(const typename hp::MappingCollection<dim, spacedim> &,
                       const unsigned int)>>) = default;

    /**
     * Move assignment operator.
     */
    MappingCollection<dim, spacedim> &
    operator=(MappingCollection<dim, spacedim> &&) = default; // NOLINT

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
  };


  /**
   * Many places in the library by default use (bi-,tri-)linear mappings
   * unless users explicitly provide a different mapping to use. In these
   * cases, the called function has to create a $Q_1$ mapping object, i.e., an
   * object of kind MappingQ(1). This is costly. It would also be
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
  template <class... MappingTypes>
  MappingCollection<dim, spacedim>::MappingCollection(
    const MappingTypes &...mappings)
  {
    static_assert(
      is_base_of_all<Mapping<dim, spacedim>, MappingTypes...>::value,
      "Not all of the input arguments of this function "
      "are derived from Mapping<dim, spacedim>!");

    // loop over all of the given arguments and add the mappings to
    // this collection. Inlining the definition of mapping_pointers causes
    // internal compiler errors on GCC 7.1.1 so we define it separately:
    const auto mapping_pointers = {
      (static_cast<const Mapping<dim, spacedim> *>(&mappings))...};
    for (const auto p : mapping_pointers)
      push_back(*p);
  }



  template <int dim, int spacedim>
  MappingCollection<dim, spacedim>
    StaticMappingQ1<dim, spacedim>::mapping_collection =
      MappingCollection<dim, spacedim>(MappingQ1<dim, spacedim>{});


#ifndef DOXYGEN
// Declare the existence of explicit instantiations of the class
// above. This is not strictly necessary, but tells the compiler to
// avoid instantiating templates that we know are instantiated in
// .cc files and so can be referenced without implicit
// instantiations. In the current case, this also avoids certain
// warnings issues by clang and newer (LLVM-based) Intel compilers.
//
// Unfortunately, this does not seem to work when building modules
// because the compiler (well, Clang at least) then just doesn't
// instantiate these classes at all, even though their members are
// defined and explicitly instantiated in a .cc file.
#  ifndef DEAL_II_BUILDING_CXX20_MODULE
  extern template struct StaticMappingQ1<1, 1>;
  extern template struct StaticMappingQ1<1, 2>;
  extern template struct StaticMappingQ1<1, 3>;
  extern template struct StaticMappingQ1<2, 2>;
  extern template struct StaticMappingQ1<2, 3>;
  extern template struct StaticMappingQ1<3, 3>;

#    ifndef _MSC_VER
  extern template MappingCollection<1, 1>
    StaticMappingQ1<1, 1>::mapping_collection;
  extern template MappingCollection<1, 2>
    StaticMappingQ1<1, 2>::mapping_collection;
  extern template MappingCollection<1, 3>
    StaticMappingQ1<1, 3>::mapping_collection;
  extern template MappingCollection<2, 2>
    StaticMappingQ1<2, 2>::mapping_collection;
  extern template MappingCollection<2, 3>
    StaticMappingQ1<2, 3>::mapping_collection;
  extern template MappingCollection<3, 3>
    StaticMappingQ1<3, 3>::mapping_collection;
#    endif
#  endif
#endif

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
